/** glp_oracle.c  --  vertex separation oracle **/

/***********************************************************************
 * This code is part of INNER, a linear multiobjective problem solver.
 *
 * Copyright (2016) Laszlo Csirmaz, Central European University, Budapest
 *
 * This program is free, open-source software. You may redistribute it
 * and/or modify under the terms of the GNU General Public License (GPL).
 *
 * There is ABSOLUTELY NO WARRANTY, use at your own risk.
 ***********************************************************************/
         
/*
* This part of INNER realizes the vertex separation oracle using a 
* "patched" version of glpk (Gnu Linear Program Kit), which allows to 
* specify an "ordering cone", and the solution to the LP returned is 
* the minimal one with respect to this cone.
* The VertexOracleData structure is used to communicate with the
* calling routine.
*/

#include "glp_oracle.h"
#include "report.h"
#include "params.h"

VertexOracle_t VertexOracleData;

/**********************************************************************
* The problem dimensions, question and answer space
*   int vcols   number of columns in the constraint matrix
*   int vrows   number of rows in the constraint matrix
*   int vobjs   dimension of the objective space: number of objectives
*   double vfacet[vobjs]
*               the question direction
*   double vvertex[vobjs]
*               the answer vertex which minimizes the facet direction
*/

#define vcols	PARAMS(ProblemColumns)
#define vrows	PARAMS(ProblemRows)
#define vobjs	PARAMS(ProblemObjects)
#define vfacet	VertexOracleData.ofacet
#define vvertex VertexOracleData.overtex

/*==================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CSL /* define CSL to allow the modifications */
#include "glpk.h"

/**********************************************************************
* The glpk LP object and parameter structure
*   glp_prob *P    the glpk problem object
*   glp_smcp parm  the glpk parameter structure
*/

static glp_prob *P=NULL;
static glp_smcp parm; /* glp parameter structure */

/**********************************************************************
* Storage for the constraint matrix, the objectives and the shuffle
*   arrays. The matrix and the shuffle arrays will be freed after
*   loading the problem into 'P'.
*
* double M(row,col)    temporary storage for the constraint matrix;
*                      indices go from 1
* double OBJ(obj,col)  the objective matrix; indices go from 1
* int vlp_rowidx[1:rows] shuffling rows, temporary
* int vlp_colidx[1:cols] shuffling columns, temporary
* int vlp_objidx[1:objs] shuffling objectives, temporary
*
* int allocate_vlp(row,cols,objs)
*   allocate memory for vfacet, vvertex, OBJ, and temporal storage.
*   Return value: -1: out of memory (should not happen)
*/

static double *vlp_M;		/* temporary storage for M */
static double *vlp_OBJ;		/* the objectives, keep it here */
static int *vlp_rowidx;		/* temporary row indices */
static int *vlp_colidx;		/* column permutation */
static int *vlp_objidx;		/* object permutation */

#define M(r,c)		vlp_M[(r)+((c)-1)*vrows-1]
#define OBJ(o,c)	vlp_OBJ[(o)+((c)-1)*vobjs-1]

#define xalloc(ptr,type,size)	\
    (ptr=(type*)calloc(size,sizeof(type)))==NULL

static int allocate_vlp(int rows, int cols, int objs)
{   vrows=rows; vcols=cols; vobjs=objs;
    if(xalloc(vfacet,double,objs+1) ||
       xalloc(vvertex,double,objs+1) ||
       xalloc(vlp_OBJ,double,objs*cols) ||
       xalloc(vlp_M,double,rows*cols) ||
       xalloc(vlp_rowidx,int,rows+1) ||
       xalloc(vlp_colidx,int,cols+1) ||
       xalloc(vlp_objidx,int,objs+1))
         return -1; // out of memory
    return 0;
}

/**********************************************************************
* Make a random permutation of an array
*
* void perm_array(int len, int array[1..len])
*   make a random permutation of the array
*/
static inline int mrandom(int v)
{ return v<=1?0 : (random()&0x3fffffff)%v; }
static inline void perm_array(int len, int arr[/* 1:len */])
{int i,j,t;
    for(i=1;i<len-1;i++){
        j=i+mrandom(len-i);
        t=arr[i];arr[i]=arr[j];arr[j]=t;
    }
}

/**********************************************************************
* Read the constraint matrix and objectives from a vlp file
*
* int MAX_LINELEN = 80
*   maximum line length expected in a vlp file
* char inpline[MAX_LINELEN]
*   the next line read form the vlp file
*
* int nextline(FILE *f)
*   read the next line to inpline[] from the open stream 'f'. Ignore
*       leading spaces and merge spaces. Convert A-Z to a-z.
*   Return value:
*        1: next line is in inpline[]
*        0: EOF
*
* int vlp_type_ok(char ctrl,int parno)
*   checks the ctrl char against the supplied number of args:
*    'f' (free) none
*    'u','l','s' (upper,lower, fixed) one
*    'd' (both) two
*
* int glp_type(char ctrl)
*   returns the glpk version of the ctrl char, or -1 if error.
*
* int read_vlp()
*   open and read the problem from the vlp file. Return value:
*     0: file read, memory allocated, dimensions stored,
*        glpk LP object P is initialized and ready to use.
*     1: some error; errors are reported as R_fatal. The vlp file
*        is not necessarily closed; the program should abort.
*/

/*------------------------------------------------------------------*/
/* read a single line from a VLP file */
#define MAX_LINELEN	80 /* maximal length of a file */

static char inpline[MAX_LINELEN+1]; /* contains the next vlp line */

/** read the next vlp line to inpline **/
static int nextline(FILE *f)
{int i,sp,ch;
    i=0;sp=0; memset(inpline,0,MAX_LINELEN+1);
    while((ch=getc(f))>=0){
        if(ch=='\n'){ if(i==0){sp=0; continue;} return 1; }
        if(ch==' '||ch=='\t'){sp=1; continue;}
        if(ch<=0x20 || ch>126) continue; /* ignore these characters */
        if(sp && i>0){if (i<MAX_LINELEN){inpline[i]=' '; i++;} }
        if('A'<=ch && ch<='Z') ch +='a'-'A'; /* upper case => lower case */
        sp=0; if(i<MAX_LINELEN){inpline[i]=ch; i++; }
    }
    /** EOF **/
    return i>0?1:0;
}

static int vlp_type_ok(char ctrl, int parno)
{   switch(ctrl){
  case 'f': return parno==2;
  case 'u': case 'l': case 's': return parno==3;
  case 'd': return parno==4;
    }
    return 0;
}
static int glp_type(char ctrl)
{   switch(ctrl){
  case 'f': return GLP_FR;
  case 'u': return GLP_UP;
  case 'l': return GLP_LO;
  case 's': return GLP_FX;
  case 'd': return GLP_DB;
    }
    return -1;
}
/* read a vlp problem from a file */
int read_vlp(void)
{FILE *f; int rows,cols,objs; int i,j,cnt; double p,b1,b2; char ctrl;
 double dir=1.0;
    f=fopen(PARAMS(VlpFile),"r");
    if(!f){
        report(R_fatal,"Cannot open vlp file %s for reading\n",PARAMS(VlpFile));
        return 1;
    }
    rows=cols=0;
    P = glp_create_prob();
    while(nextline(f)) switch(inpline[0]){ /* read next input line */
       case 'c': // comment line, print out nonempty comment lines before the first p line
                 if(rows==0 && inpline[1]){
                     report(R_warn,"C%s\n",inpline+1); 
                 }
                 continue;
       case 'e': break; // end
       case 'p': if(rows>0){
                    report(R_fatal,"read_vlp: second p line in %s:\n   %s\n",
                              PARAMS(VlpFile),inpline); return 1; 
                 }
                 dir=+1.0; PARAMS(Direction)=0;
                                    /**   <rows> <cols> <> <objs> <> **/
                 cnt=sscanf(inpline,"p vlp min %d %d %*d %d %*d",&rows,&cols,&objs);
                 if(cnt==0){
                    dir=-1.0; PARAMS(Direction)=1;
                    cnt=sscanf(inpline,"p vlp max %d %d %*d %d %*d",&rows,&cols,&objs);
                 }
                 if(cnt!=3 || rows<=1 || cols<=1 || objs<1 ){
                    report(R_fatal,"read_vlp: wrong p line in %s\n   %s\n",
                               PARAMS(VlpFile),inpline); return 1;
                 }
                 if(allocate_vlp(rows,cols,objs)){
                     report(R_fatal,"read_vlp: out of memory for %s:\n   %s\n",
                               PARAMS(VlpFile),inpline); return 1;
                 }
                 for(i=0;i<=rows;i++) vlp_rowidx[i]=i;
                 for(j=0;j<=cols;j++) vlp_colidx[j]=j;
                 if(PARAMS(ShuffleMatrix)){
                     perm_array(rows,vlp_rowidx);
                     perm_array(cols,vlp_colidx);
                 }
		 glp_add_cols(P,cols); glp_add_rows(P,rows);
                 glp_set_multiobj_number(P,objs);
                 continue;
       case 'j': // j <col> [ f || l <val> | u <val> | d <val1> <val2> | s <val> ]
                 if(rows==0){
                    report(R_fatal,"read_vlp: j line before p in %s\n  %s\n",
                               PARAMS(VlpFile),inpline); return  1;
                 }
                 b1=b2=0.0; 
                 cnt=sscanf(inpline,"j %d %c %lg %lg",&j,&ctrl,&b1,&b2);
                 if(cnt<2 || cols<j || cols<1 || !vlp_type_ok(ctrl,cnt)){
                    report(R_fatal,"read_vlp: wrong j line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 if(cnt<4) b2=b1; // GLP_UP uses b2 as the value
                 glp_set_col_bnds(P,vlp_colidx[j],glp_type(ctrl),b1,b2);
                 continue;
       case 'i': if(rows==0){
                    report(R_fatal,"read_vlp: i line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 // i <row> [ f | l <val> | u <val> | d <val1> <val2> | s <val> ]
                 b1=b2=0.0;
                 cnt=sscanf(inpline,"i %d %c %lg %lg",&i,&ctrl,&b1,&b2);
                 if(cnt<2 || rows<i || i<1 || !vlp_type_ok(ctrl,cnt)){
                    report(R_fatal,"read_vlp: wrong i line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 if(cnt<4) b2=b1; // GLP_UP uses b2 as the value
                 glp_set_row_bnds(P,vlp_rowidx[i],glp_type(ctrl),b1,b2);
                 continue;
       case 'a': if(rows==0){
                    report(R_fatal,"read_vlp: a line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 cnt=sscanf(inpline,"a %d %d %lg",&i,&j,&p);
                 if(cnt!=3 || rows<i || i<1 || cols<j || j<1){
                    report(R_fatal,"read_vlp: wrong a line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 M(vlp_rowidx[i],vlp_colidx[j])=p; // store it
                 continue;
       case 'o': if(rows==0){
                    report(R_fatal,"read_vlp: o line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 cnt=sscanf(inpline,"o %d %d %lg",&i,&j,&p);
                 if(cnt!=3 || objs<i || i<1 || cols<j || j<1){
                    report(R_fatal,"read_vlp: wrong o line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 OBJ(i,vlp_colidx[j])=dir*p; // store it
                 continue;
       default: report(R_fatal,"read_vlp: unknown line in %s\n  %s\n",
                           PARAMS(VlpFile),inpline); return 1;
    }
    fclose(f);
    if(rows==0){
       report(R_fatal,"read_vlp: no 'p' line in %s\n",PARAMS(VlpFile)); return 1; 
    }
    free(vlp_colidx);
    // upload constraints into P
    for(i=0;i<=rows;i++) vlp_rowidx[i]=i;
    for(j=1;j<=cols;j++){
        glp_set_mat_col(P,j,rows,vlp_rowidx,&M(1,j)-1);
    }
    free(vlp_M); free(vlp_rowidx);
    // upload objectives into P
    for(i=1;i<=objs;i++) vlp_objidx[i]=i;
    if(PARAMS(ShuffleMatrix)) perm_array(objs,vlp_objidx);
    for(i=1;i<=objs;i++){
        for(j=1;j<=cols;j++){
           glp_set_multiobj_coef(P,i,j,OBJ(vlp_objidx[i],j));
        }
    }
    free(vlp_objidx);
    glp_set_obj_dir(P,GLP_MIN); // minimize
    // suppress glpk message if not verbose ...
    if(PARAMS(OracleMessage)<2) glp_term_out(GLP_OFF);
    glp_sort_matrix(P);
    if(PARAMS(OracleScale)) glp_scale_prob(P,GLP_SF_AUTO);
    glp_adv_basis(P,0);		// make this optimization
    glp_term_out(GLP_ON);
//    set_oracle_parameters(); will be set on the main program
    return 0;
}

/**********************************************************************
* Set the LP solver parameters from the configuration.
*
* void set_oracle_parameters(void)
*   verbosity: 0: no, 1: error; 2: normal, 3: verbose
*   output frequency: indicate that the LP solver is working
*   method: primal, dual
*   pricing: standard or steepest edge
*   ratio test: standard of Harris
*   iteration limit
*   time limit (in seconds)
*/

void set_oracle_parameters(void)
{   glp_init_smcp(&parm);
    // multiobjective
    if(PARAMS(ProblemObjects)>1) parm.mobj=GLP_ON;
    // verbosity
    switch(PARAMS(OracleMessage)){
      case  0: parm.msg_lev=GLP_MSG_OFF; break;
      case  1: parm.msg_lev=GLP_MSG_ERR; break;
      case  2: parm.msg_lev=GLP_MSG_ON; break;
      default: parm.msg_lev=GLP_MSG_ALL; break;
    }
    // method, pricing, ratio test
    parm.meth = GLP_PRIMAL;	// PRIMAL,DUAL,DUALP
    if(PARAMS(OracleMethod)) parm.meth=GLP_DUAL;
    parm.pricing = GLP_PT_STD;
    if(PARAMS(OraclePricing)) parm.pricing=GLP_PT_PSE;
    parm.r_test = GLP_RT_STD;	// HAR
    if(PARAMS(OracleRatioTest)) parm.r_test=GLP_RT_HAR;
    // iteration and time limit
    parm.it_lim = 100000;	// iteration limit
    if(PARAMS(OracleItLimit)>=1000) parm.it_lim=PARAMS(OracleItLimit);
    if(PARAMS(OracleItLimit)==0) parm.it_lim=0; // no limit
    parm.tm_lim = 10000;	// time limit in milliseconds
    if(PARAMS(OracleTimeLimit)>=5) parm.tm_lim=1000*PARAMS(OracleTimeLimit);
    if(PARAMS(OracleTimeLimit)==0) parm.tm_lim=0; // no limit
}

/**********************************************************************
* Ask the oracle
*
* int ask_oracle(void)
*    the question and the answer is provided in VertexOracleData.
*
* Return values:
*   ORACLE_OK     the minimal vertex is stored in VertexOracleData.
*                 coordinates are rounded to the nearest rational value
*                 when "RoundVertices" is set.
*   ORACLE_UNBND  the polytope is not bounded from below
*   ORACLE_EMPTY  the polytope is empty (no feasible solution)
*   ORACLE_LIMIT  either time or iteration limit is reached
*   ORACLE_FAIL   the LP solver failed to solve the problem
* Errors are reported as R_err.
*
* void round_to_int(double *v)
*  v is expanded as a continued fraction using three iterations.
*  If the error is smaller than PARAMS(RoundEps), then v is
*  replaced by the value of the fraction.
* char *glp_status_msg(int code)
* char *glp_return_msg(int code)
*  return the verbatim error message corresponding to the glpk
*  code (to be printed out).
*/
#define ROUND_EPS	PARAMS(RoundEps)

static inline long intfloor(double d)
{ return d<0.0 ? (long)(d-0.5) : (long)(d+0.5); 
}

static inline void round_to_int(double *v)
{double ip,iip,v2,ip2,ip3;
    ip=intfloor(*v); iip=*v-ip;
//*v  (-2.5,-1.5]  (-1.5,-0.5] (-0.5,0] [0,0.5) [0.5,1.5)  [1.5,2.5)
//ip    -2          -1           0          0      1          2
//iip (-0.5,0.5]   (-0.5,0.5]  (-0.5,0] [0.0.5] [-0.5,0.5) [-0.5,0.5)
    if(-ROUND_EPS<iip && iip<ROUND_EPS){ *v=ip; return; }
    // *v==ip+iip, -0.50<=iip<=0.5
    v2=1.0/iip; ip2=intfloor(v2); iip=v2-ip2;
    if(-ROUND_EPS<iip && iip<ROUND_EPS){ *v=ip+1.0/ip2; return; }
    v2=1.0/iip; ip3=intfloor(v2); iip=v2-ip3;
    if(-ROUND_EPS<iip && iip<ROUND_EPS){ *v=ip+1.0/(ip2+1.0/ip3); return; }
}

static char *glp_status_msg(int stat)
{static char *statmsg[] = {
"the problem is undefined",		// GLP_UNDEF
"solution is feasible",			// GLP_FEAS
"solution is infeasible",		// GLP_INFEAS
"the problem has no feasible solution",	// GLP_NOFEAS
"solution is optimal",			// GLP_OPT
"the problem is unbounded",		// GLP_UNBND
};
    if(1<=stat && stat<=6) return statmsg[stat-1];
    return "unknown solution status";
}

static char *glp_return_msg(int retval)
{static char *retmsg[] = {
"invalid basis",			// GLP_EBADB	 *
"singular matrix",			// GLP_ESING	 *
"ill-conditioned matrix",		// GLP_ECOND	 *
"invalid bounds",			// GLP_EBOUND
"solver failed",			// GLP_EFAIL	 *
"objective lower limit reached",	// GLP_EOBJLL
"objective upper limit reached",	// GLP_EOBJUL
"iteration limit exceeded",		// GLP_EITLIM	 *
"time limit exceeded",			// GLP_ETMLIM	 *
"no primal feasible solution",		// GLP_ENOPFS
"no dual feasible solution",		// GLP_ENODFS
"root LP optimum not provided",		// GLP_EROOT
"search terminated by application",	// GLP_ESTOP
"relative mip gap tolerance reached",	// GLP_EMIPGAP
"no primal/dual feasible solution",	// GLP_ENOFEAS
"no convergence",			// GLP_ENOCVG
"numerical instability",		// GLP_EINSTAB
"invalid data",				// GLP_EDATA
"result out of range",			// GLP_ERANGE
};
    if(1<=retval && retval<=0x13) return retmsg[retval-1];
    return "unknown error";
}

/** the question is in VertexOracleData.ofacet;
    the answer goes to VertexOracleData.overtex
    return value: 
      ORACLE_OK    : OK
      ORACLE_LIMIT : resource limit reached
      ORACLE_FAIL  : some problem,
**/
int ask_oracle(void)
{int i,j,ret; double d;
    /* set up the problem to be minimized */
    for(j=1;j<=vcols;j++){
        d=0.0;
        for(i=1;i<=vobjs;i++) d += OBJ(i,j)*vfacet[i-1];
        glp_set_obj_coef(P,j,d);
    }
    ret=glp_simplex(P,&parm);
    if(ret==GLP_EFAIL){ /* give it another change */
        if(PARAMS(OracleMessage)<2) glp_term_out(GLP_OFF);
        glp_adv_basis(P,0);
        glp_term_out(GLP_ON);
        ret=glp_simplex(P,&parm);
    }
    if(ret){
        report(R_fatal,"The oracle says: %s (%d)\n",glp_return_msg(ret),ret);
        // one can continue if  ret==GLP_EITLIM || ret==GLP_ETMLIM
        return (ret==GLP_EITLIM || ret==GLP_ETMLIM)? ORACLE_LIMIT : ORACLE_FAIL; 
    }
    ret=glp_get_status(P);
    if(ret != GLP_OPT){
        report(R_fatal,"The oracle says: %s (%d)\n",glp_status_msg(ret),ret);
        return ret==GLP_NOFEAS? ORACLE_EMPTY : 
               ret==GLP_UNBND ? ORACLE_UNBND : ORACLE_FAIL;
    }
    /* recover the solution vertex */
    for(i=1;i<=vobjs;i++){
        d=0.0;
        for(j=1;j<=vcols;j++) d+=OBJ(i,j)*glp_get_col_prim(P,j);
        if(PARAMS(RoundVertices))round_to_int(&d);
        vvertex[i-1]=d;
    }
    return ORACLE_OK;
}

/**********************************************************************
* Get oracle iteration count
*
* int get_oracle_rounds(void)
*   use the undocumented glpk call to get the number of iterations.
*   Used in statistics.
*/

int get_oracle_rounds(void)
{   return glp_get_it_cnt(P); }

/* EOF */

