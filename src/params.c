/** params.c -- global parameters **/

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

#include "params.h"
#include "report.h"
#include "version.h"
#include <stdio.h>

/***********************************************************************
* Default values in the config file
*/
/* oracle parameters */
#define DEF_OracleMessage	1	/* error */
#define DEF_OracleItLimit	10000
#define DEF_OracleTimeLimit	20	/* in seconds */
#define DEF_OracleMethod	0	/* PRIMAL/DUAL */
#define DEF_OraclePricing	1	/* STD / steepest */
#define DEF_OracleRatioTest	1	/* STD / Harris */
#define DEF_OracleScale		1	/* scale */
#define DEF_ShuffleMatrix	1
#define DEF_RoundVertices	1
/* DD parameters */
#define DEF_RandomFacet		1	/* yes */
#define DEF_ExactFacetEq	0	/* no */
#define DEF_RecalculateFacets	100
#define DEF_CheckConsistency	0
#define DEF_ExtractAfterBreak	1	/* yes */
/* vertex pool */
#define DEF_VertexPoolSize	0	/* don't use vertex pool */
#define DEF_OracleCallLimit	1	/* stop after the first unsuccessful call */
/* number of threads */
#define DEF_Threads		0	/* number of threads */
/* randomness */
#define DEF_TrueRandom		1	/* yes */
/* Tolerances */
#define DEF_RoundEps		1e-9
#define DEF_ScaleEps		3e-9
#define DEF_PolytopeEps		1.3e-8
#define DEF_LineqEps		8e-8	/* 6.0*PolytopeEps */
#define DEF_FacetRecalcEps	1e-6
/* Reporting */
#define DEF_MessageLevel	2	/* on */
#define DEF_PrintParams		1	/* yes */
#define DEF_PrintStatistics	1	/* yes */
#define DEF_ProgressReport	5	/* in seconds */
#define DEF_MemoryReport	0	/* don't report */
#define DEF_VertexReport	1	/* yes */
#define DEF_VertexAsFraction	1	/* yes */
#define DEF_PrintVertices	2	/* partial results */
#define DEF_PrintFacets		0	/* don't report */
#define DEF_SaveVertices	2	/* partial */
#define DEF_SaveFacets		1	/* on normal exit only */
/* name of this program */
#ifndef PROG
  #define PROG			inner
#endif
#define PROGNAME		mkstringof(PROG)

struct params_t GlobalParams;

/***********************************************************************
* Content of the default config file
* void dump_config(void)
*    dump a config file with the default values
*/

#define CFG(x,y)   " " #x " = " mkstringof(DEF_##x) "   # " y "\n"

#define BOOL	"0 = no, 1 = yes"
#define BOOL2	"0 = no, 1 = on normal exit only, 2 = always"
#define INTEGER	"non-negative integer"
#define POSINT	"positive integer"
#define REAL	"positive real number"


#define DEFAULT_CONFIG_FILE	\
"# " PROGNAME ".cfg -- configuration file for the program " PROGNAME "\n"\
"\n"\
"# This file defines different parameters of the MOLP solver. Values\n"\
"# defined here are superseded by command line options. Unspecified\n"\
"# keywords take their default values.\n"\
"# In this file everything after a # symbol is ignored.\n"\
"\n"\
"##########################\n"\
"#  ALGORITHM parameters  #\n"\
"##########################\n"\
"#\n"\
CFG( RandomFacet, BOOL)  \
"#    pick the next facet to be asked the oracle about randomly.\n"\
"#\n"\
CFG( ExactFacetEq, BOOL) \
"#    when a facet is created, recompute its equation immediately\n"\
"#    from the set of adjacent vertices.\n"\
"#\n"\
CFG( RecalculateFacets, INTEGER) \
"#    after that many iterations recalculate all facet equations\n"\
"#    from the set of adjacent vertices. Should be zero (never),\n"\
"#    or at least 5.\n"\
"#\n"\
CFG( CheckConsistency, INTEGER) \
"#    after that many iterations check the consistency of the data\n"\
"#    structure against numerical errors. Should be zero (never),\n"\
"#    or at least 5.\n"\
"#\n"\
CFG( ExtractAfterBreak, BOOL) \
"#    when the program receives a " mkstringof(INNER_SIGNAL) " signal, continue extracting\n"\
"#    new vertices by asking the oracle about every facet of the\n"\
"#    actual approximating polyhedron. Can be very time consuming.\n"\
"#    Second signal aborts post-processing.\n"\
"#\n"\
CFG( VertexPoolSize, INTEGER) \
"#    size of the vertex pool; add the vertex to the approximation\n"\
"#    which discards the largest number of existing facets. Should\n"\
"#    be zero (don't use it) or at least 5. Using vertex pool makes\n"\
"#    additional oracle calls, but can simplify the approximating\n"\
"#    polytopes. See OracleCallLimit. Should be less than " mkstringof(MAX_VERTEX_POOL) ".\n"\
"#\n"\
CFG( OracleCallLimit, INTEGER ) \
"#    the maximal number of unsuccessful oracle calls during an\n"\
"#    iteration when filling the vertex pool. Zero means no limit;\n"\
"#    otherwise should be less than " mkstringof(MAX_OCALL_LIMIT)  ".\n"\
"#\n"\
CFG( Threads, INTEGER) \
"#    number of threads to use; should be less than " mkstringof(MAX_THREADS) ". Zero means\n"\
"#    use as many as are available; 1 means don't use threads.\n"\
"#\n"\
"##########################\n"\
"#   ORACLE parameters    #\n"\
"##########################\n"\
"#\n"\
CFG( OracleMessage, "0 = quiet, 1 = error, 2 = on, 3 = verbose") \
"#    oracle (glpk) message level.\n"\
"#\n"\
CFG( OracleMethod, "0 = primal, 1 = dual") \
"#    the LP method used by the oracle.\n"\
"#\n"\
CFG( OraclePricing, "0 = standard, 1 = steepest edge") \
"#    the LP pricing method.\n"\
"#\n"\
CFG( OracleRatioTest, "0 = standard, 1 = Harris' two pass") \
"#    the LP ratio test.\n"\
"#\n"\
CFG( OracleTimeLimit, INTEGER) \
"#    time limit for each oracle call in seconds, 0 = unlimited.\n"\
"#\n"\
CFG( OracleItLimit, INTEGER) \
"#    iteration limit for each oracle call, 0 = unlimited.\n"\
"#\n"\
CFG( OracleScale, BOOL) \
"#    scale the constraint matrix; helps numerical stability.\n"\
"#\n"\
CFG( ShuffleMatrix, BOOL) \
"#    shuffle the rows and columns of the constraint matrix randomly.\n"\
"#\n"\
CFG( RoundVertices, BOOL) \
"#    when the oracle reports a result vertex, round its coordinates\n"\
"#    to the nearest rational with small denominator.\n"\
"#\n"\
"##########################\n"\
"#       REPORTING        #\n"\
"##########################\n"\
"#\n"\
CFG( MessageLevel, "0 = quiet, 1 = error, 2 = all, 3 = verbose")\
"#    report level, quiet means no messages at all. Command line\n"\
"#    option -m[0..3] overrides this value.\n"\
"#\n"\
CFG( ProgressReport, INTEGER) \
"#    minimum time between two progress reports (in seconds). Should\n"\
"#    be zero (no progress report), or at least 5. Use command line\n"\
"#    option -p <seconds> to override this value.\n"\
"#\n"\
CFG( VertexReport, BOOL) \
"#    print out each vertex (extremal solution) immediately as it is\n"\
"#    found. Use command line option -y+ (yes) or -y- (no) to override\n"\
"#    the value defined here.\n"\
"#\n"\
CFG( MemoryReport, BOOL) \
"#    report the size and location, whenever it changes, of memory\n"\
"#    blocks storing the data structure.\n"\
"#\n"\
CFG( VertexAsFraction, BOOL) \
"#    if possible, print (and save) vertex coordinates as fractions\n"\
"#    rather than floating point numerals.\n"\
"#\n"\
CFG( PrintStatistics, BOOL) \
"#    print out resources used (number of iterations, ridge tests, etc.)\n"\
"#\n"\
CFG( PrintParams, BOOL) \
"#    print out parameter values which are not the default ones.\n"\
"#\n"\
CFG( PrintVertices, BOOL2) \
"#    print out all known vertices when the program terminates.\n"\
"#\n"\
CFG( PrintFacets, BOOL2) \
"#    print out all known facets when the program terminates.\n"\
"#\n"\
CFG( SaveVertices, BOOL2) \
"#    when the program terminates save known vertices to the file\n"\
"#    specified after the command line option '-o'. For the file\n"\
"#    specified after '-ov' both 0 and 1 means \"save on normal\n"\
"#    exit only\".\n"\
"#\n"\
CFG( SaveFacets, BOOL2) \
"#    when the program terminates save known facets to the file\n"\
"#    specified after command line options '-o'. For the file\n"\
"#    specified after '-of' both 0 and 1 means \"save on normal\n"\
"#    exit only\".\n"\
"#\n"\
"##########################\n"\
"#       TOLERANCES       #\n"\
"##########################\n"\
"#\n"\
"#  >>> Change these values with great care ...      <<<\n"\
"#  >>> ... and don't forget to delete the leading # <<<\n"\
"#\n"\
"#" CFG( PolytopeEps, REAL) \
"#    a facet and a vertex are considered adjacent if their distance\n"\
"#    is smaller than this value.\n"\
"#\n"\
"#" CFG( ScaleEps, REAL) \
"#    coefficients in the scaled facet equation are rounded to the\n"\
"#    nearest integer if they are closer to it than this value.\n"\
"#\n"\
"#" CFG( LineqEps, REAL) \
"#    when solving a system of linear equations for a facet equation,\n"\
"#    a coefficient smaller than this is considered to be zero.\n"\
"#\n"\
"#" CFG( RoundEps, REAL) \
"#    if Oracle vertex reports are rounded (RoundVertices=1), the\n"\
"#    tolerance in the rounding algorithm.\n"\
"#\n"\
"#" CFG( FacetRecalcEps, REAL) \
"#    when recalculating facet equations, report numerical instability\n"\
"#    if the new and old coordinates differ at least that much.\n"\
"#\n"\
"# *** end of " PROGNAME ".cfg ***\n\n"

static void dump_config(void)
{     printf("%s", DEFAULT_CONFIG_FILE); }

#undef CFG
#undef BOOL
#undef BOOL2
#undef INTEGER
#undef POSINT
#undef REAL
#undef DEFAULT_CONFIG_FILE

/***********************************************************************
* Help routines
* void short_help(void)
*    show a short help
* void long_help(void)
*    show all command line options
* void vlp_help(void)
*    description of the vlp file
* void out_help(void)
*    how vertices and facets are printed 
* void version(void)
*    version and copyright
*/

static void short_help(void) {printf(
"Solve a multiobjective linear program using the inner approximation method.\n"
"Usage: " PROGNAME " [options] <vlp file>\n"
"Some of the options are:\n"
"  -h               display this short help\n"
"  --help           display all options\n"
"  -c <config-file> specify configuration file\n"
"  -o <file>        save the solution to <file>\n"
"  -q               quiet, no messages, no statistics\n"
"  -p0              no progress report\n"
"  -y-              do not report vertices when generated\n"
"Previous content of the output file is deleted without warning.\n"
COPYRIGHT "\n"
);}

static void long_help(void){ printf(
"Solve a multiobjective linear program using the inner approximation method.\n"
"Usage: " PROGNAME " [options] <vlp file>\n"
"Options are:\n"
"  -h               display a short help\n"
"  --help           display all options\n"
"  --help=vlp       describe vlp file format\n"
"  --help=out       describe output format\n"
"  --version        version and copyright information\n"
"  --dump           dump the default config file and quit\n"
"  --config=<config-file>\n"
"  -c <config-file> read configuration from the given file (see --dump)\n"
"  -o <file>        save result (both vertices and facets) to <file>\n"
"  -ov <file>       save vertices to <file>\n"
"  -of <file>       save facets to <file>\n"
"  --name=NAME\n"
"  -n NAME          specify the problem name\n"
"  -m[0..3]         set message level: 0: none, 3: verbose\n"
"  -q               quiet, same as -m0. Implies --PrintStatistics=0\n"
"  -p T             progress report in every T seconds (default: T=5)\n"
"  -p 0             no progress report\n"
"  -y+              report vertices immediately when generated (default)\n"
"  -y-              do not report vertices when generated\n"
"  --KEYWORD=value  change value of a config keyword (see --dump)\n"
"Previous content of output files are deleted without warning.\n"
COPYRIGHT "\n"
);}

static void vlp_help(void) {printf(
"****************************\n"
"***   VLP input format   ***\n"
"****************************\n"
"\n"
"A multiple objective linear program (MOLP) is a linear program with\n"
"multiple objective functions. The VLP format describes such a problem\n"
"in a plain text file. In this file each line begins with a lower-case\n"
"letter which identifies the line type, which can be one of\n"
"    c    comment\n"
"    p    program line, it should be the first non-comment line\n"
"    i    constraint matrix row descriptor\n"
"    j    constraint matrix column descriptor\n"
"    a    constraint matrix coefficient\n"
"    o    objective coefficient\n"
"    e    end of data, last processed line in the vlp file\n"
"Comment lines are ignored. The 'p' program line has the format\n"
"    p vlp DIR ROWS COLS ALINES OBJS OLINES\n"
"DIR is either 'min' or 'max' defining whether the problem is to minimize\n"
"or maximize the objectives. Other fields are positive integers: ROWS, COLS\n"
"are the number of rows and columns of the constraint matrix; OBJS is the\n"
"number of objectives. ALINES and OLINES are the number of 'a' and 'o'\n"
"lines in the vlp file; these numbers are ignored by this program.\n"
"\n"
"PLEASE NOTE: rows, columns and objectives ARE INDEXED STARTING FROM 1.\n"
"\n"
"A row descriptor line starting with 'i' can be one of the following:\n"
"    i ROW f            row is free, there is no constraint\n"
"    i ROW l VAL        row with lower bound, row's value is >= VAL\n"
"    i ROW u VAL        row with upper bound, row's value is <= VAL\n"
"    i ROW d VAL1 VAL2  doubly bounded row: VAL1 <= row's value <= VAL2\n"
"    i ROW s VAL        row's value is fixed to be equal to VAL\n"
"here ROW is the row's index (between 1 and the number of rows), and VAL is\n"
"a floating point constant. The form of a 'j' line is similar, it describes\n"
"the variable types: free, with lower and/or upper bound, or fixed. The\n"
"default row type is 'free', the default column type is 'fixed' with value\n"
"zero.\n"
"Elements of the constraint matrix are specified as\n"
"    a ROW COL VAL      both ROW and COL are positive integers\n"
"while the coefficients in the OBJ-th objective function are given as\n"
"    o OBJ COL VAL\n"
"VAL is a floating point constant; 'a' and 'o' lines with zero values can\n"
"be omitted.\n"
);}

static void out_help(void) {printf(
"*************************\n"
"***   Output format   ***\n"
"*************************\n"
"\n"
"The solution of a MOLP with d objectives is a list of d-dimensional points:\n"
"the EXTREMAL solutions. These points are the vertices of the d-dimensional\n"
"polyhedron which is formed by the achievable solutions and their superset\n"
"(when the problem is minimize), or subset (when the problem is maximize).\n"
"\n"
"Extremal solutions are printed in separate lines starting with 'V' followed\n"
"by the value of the d objectives separated by spaces:\n"
"    V 0 5/2 3/4 7.123456789 -1/2\n"
"Numbers are printed as fractions with small denominator whenever possible.\n"
"To print them as floating point numerals use the '--VertexAsFraction=0'\n"
"command line option, or change this value in the default config file.\n"
"\n"
"When requested, facets of the extremal polyhedron are printed in separate\n"
"lines starting with 'F' followed by d+1 floating point numerals separated by\n"
"spaces. The facet has equation f[1]*x1+...+f[d]*xd + f[d+1]=0.\n"
"    F 13 7 7 1 0 -10\n"
"\n"
"Other non-empty lines in the output start with C, and contain information\n"
"such as the name and size of the problem; whether it is a partial list; and\n"
"the number of vertices and facets generated.\n"
);}

#include "glpk.h"

static void version(void) {printf(
"This is '" PROGNAME "' Version " VERSION_STRING ", a multiobjective linear program solver,\n"
"using a patched version of glpk " mkstringof( GLP_MAJOR_VERSION.GLP_MINOR_VERSION)
" (GNU Linear Programming Kit).\n"
COPYRIGHT "\n"
);}

/***********************************************************************
* Parameter configuration
*
* int config_error
*    set if there was an error while handling configuration
*
* struct char_params, int_params, double_params
*    information about each configuration keyword:
*       a format string to parse the keyword
*       pointer to the real configuration value
*       default value
*       maximal value (not for double keywords)
*       whether the keyword has been set previously
*/
static int config_error=0; // set if there was a config error

#include "string.h"

#define CFG(x,y)   { " " #x " = %d %c", &PARAMS(x), DEF_##x , y, 0 }

static struct char_params { 
    const char *format;    // keyword = %d
    char       *ptr;       // pointer to PARAMS()
    const char deft;       // default value when none is given
    const char max;        // maximal accepted value, 0..1, 0..2, 0..3 
    char       filled;     // filled
} CHAR_PARAMS[] = {
  CFG(MessageLevel,3),
  CFG(PrintParams,1),
  CFG(PrintStatistics,1),
  CFG(VertexAsFraction,1),
  CFG(VertexReport,1),
  CFG(MemoryReport,1),
  CFG(PrintVertices,2),
  CFG(PrintFacets,2),
  CFG(SaveVertices,2),
  CFG(SaveFacets,2),
  CFG(RandomFacet,1),
  CFG(ExactFacetEq,1),
  CFG(ExtractAfterBreak,1),
  CFG(TrueRandom,1),
  CFG(ShuffleMatrix,1),
  CFG(RoundVertices,1),
  CFG(OracleMessage,3),
  CFG(OracleScale,1),
  CFG(OracleMethod,1),
  CFG(OracleRatioTest,1),
  CFG(OraclePricing,1),
  {NULL,NULL,0,0,0}
};

#undef CFG
#define CFG(x,y)   { " " #x " = %d %c", &PARAMS(x), DEF_##x , y, 0 }

static struct int_params {
    const char *format;    // keyword = %d
    int        *ptr;       // pointer to PARAMS()
    const int  deft;       // default value
    const int  max;        // maximal value
    int        filled;     // filled
} INT_PARAMS[] = {
  CFG(ProgressReport,1000000),
  CFG(RecalculateFacets,1000000),
  CFG(CheckConsistency,1000000),
  CFG(VertexPoolSize,MAX_VERTEX_POOL),
  CFG(Threads,MAX_THREADS),
  CFG(OracleCallLimit,MAX_OCALL_LIMIT),
  CFG(OracleItLimit,10000000),
  CFG(OracleTimeLimit,1000000),
  {NULL,NULL,0,0,0}
};

#undef CFG
#define CFG(x)   { " " #x " = %lg %c", &PARAMS(x), DEF_##x , 0 }

static struct double_params {
    const char   *format;   // keyword = %lg
    double       *ptr;      // pointer to PARAMS()
    const double deft;      // default value
    int          filled;    // filled
} DOUBLE_PARAMS[] = {
  CFG(RoundEps),
  CFG(ScaleEps),
  CFG(PolytopeEps),
  CFG(LineqEps),
  CFG(FacetRecalcEps),
  {NULL,NULL,0,0}
};
#undef CFG

/************************************************************************************
* int treat_keyword(char *line)
*    check if line is of the form KEYWORD=value. If yes, and KEYWORD
*    is not filled yet, change KEYWORD's value.
*    return value:
*        0: KEYWORD not found
*        1: KEYWORD found, value OK
*       -1: KEYWORD found, value is out of range
*
* void set_default_values(void)
*    set default values to unset PARAMS.
*/
static int treat_keyword(const char *line)
{{struct char_params *p; int val; char c;// bool params
    for(p=&CHAR_PARAMS[0];p->format;p++){
        if(sscanf(line,p->format,&val,&c)==1){
            if(val<0 || val>p->max) return -1;
            if(!p->filled){ p->filled=1; *(p->ptr)=(char)val; }
            return 1;
        }
    }
 }{struct int_params *p; int val; char c;// integer params
    for(p=&INT_PARAMS[0];p->format;p++){
        if(sscanf(line,p->format,&val,&c)==1){
            if(val<0 || val > p->max) return -1;
            if(!p->filled){ p->filled=1; *(p->ptr)=val; }
            return 1;
        }
    }
 }{struct double_params *p; double val; char c; // double params
    for(p=&DOUBLE_PARAMS[0];p->format;p++){
        if(sscanf(line,p->format,&val,&c)==1){
            if(val<1.01e-15 || val>0.99) return -1;
            if(!p->filled){ p->filled=1; *(p->ptr)=val; }
            return 1;
        }
    }
 }
 return 0; // not found
}

/** set_default_values(): set default values for unset PARAMS **/
static void set_default_values(void)
{{struct char_params *p; // character/bool params
    for(p=&CHAR_PARAMS[0];p->format;p++)if(!p->filled){
        *(p->ptr)=p->deft;
    }
 }{struct int_params *p; // integer params
    for(p=&INT_PARAMS[0];p->format;p++)if(!p->filled){
        *(p->ptr)=p->deft;
    }
 }{struct double_params *p; // double params
    for(p=&DOUBLE_PARAMS[0];p->format;p++)if(!p->filled){
        *(p->ptr)=p->deft;
    }
 }
}

/**********************************************************************
* Read config file and parse its content
*
* int MAX_LINELEN
*    maximum line length read from the config file
* char inpline[MAX_LINELEN+1]
*    the next non-empty, non comment line of the config
* int nextline(FILE *f)
*    read the next line to inpline[]. Ignore leading spaces,
*      merge spaces, take care of # comments
*    Return value:
*      1: next non-empty line is in inpline[]
*      0: EOF
* void read_config_file(void)
*    if ConfigFile is specified, open and read values. Issue
*    fatal messages for unrecognized lines.
*/
#define MAX_LINELEN	80 /* maximal considered line length */
static char inpline[MAX_LINELEN+1];/* read next line here */

static int nextline(FILE *f)
{int i,sp,ch,cmt;
    i=0;sp=0;cmt=0; memset(inpline,0,MAX_LINELEN+1);
    while((ch=getc(f))>=0){
        if(ch=='\n'){ if(i==0){sp=0; cmt=0; continue;} return 1; }
        if(cmt) continue;
        if(ch==' '||ch=='\t'){sp=1; continue;}
        if(ch<=0x20 || ch>126) continue; /* ignore these characters */
        if(ch=='#'){ cmt=1; continue; }
        if(sp && i>0){if (i<MAX_LINELEN){inpline[i]=' '; i++;} }
        sp=0; if(i<MAX_LINELEN){inpline[i]=ch; i++; }
    }
    /** EOF **/
    return i>0?1:0;
}

static void read_config_file(void)
{FILE *f;
    if(PARAMS(ConfigFile)==NULL) return; // no config file
    if(!*PARAMS(ConfigFile)) return; // empty string
    f=fopen(PARAMS(ConfigFile),"r");
    if(!f){
        report(R_fatal,"Cannot open config file %s\n",PARAMS(ConfigFile));
        config_error++; return;
    }
    while(nextline(f)) switch(treat_keyword(inpline)){
  case -1: report(R_fatal,"In config file \"%s\" keyword value is out of range:\n  %s\n",
               PARAMS(ConfigFile),inpline);
           config_error++;
           break;
  case  0: report(R_fatal,"In config file \"%s\" misspelled keyword or wrong syntax:\n  %s\n",
               PARAMS(ConfigFile),inpline);
           config_error++;
           break;
  default: break; // OK
    }
    fclose(f);
}

/***********************************************************************
* Handle commandline options
*
* int handle_options(int argc, char *argv[])
*    go over all options; handle --help, --version, --dump
*    at the end read config file (if specified), and then
*    set default values for unset PARAMS().
*    Return value:
*      0:  options are OK, can continue
*      1:  option handled, exit normally
*     -1:  some error, error message issued, exit with error
*
* int integer_arg(char *a1, char *a2, int *c, int *val)
*     the option a1 needs an integer parameter. It starts at
*     the second char of a1, or, if it is \0, then in a2. In the
*     latter case increase c by one. The value is returned in
*     *val. Return value: 0: OK, 1: some error, message reported.
*/
static int integer_arg(const char *a1, const char *a2, int *c, int *val)
{int i;
    if(a1[2]==0){
        if(a2){
            if('0'<= *a2 && *a2<='9') (*c)++;
            for(i=0; '0'<=*a2 && *a2<='9'; a2++){i=i*10+(*a2)-'0';}
            if(*a2==0){ *val=i; return 0; }
        }
        report(R_fatal,"Option %s requires an integer argument.\n",a1);
        config_error++; return 1;
    } else {
        a2=a1+2;
        for(i=0; '0'<=*a2 && *a2<='9'; a2++){i=i*10+(*a2)-'0';}
        if(*a2==0){ *val=i; return 0; }
    } // followed by non-digit sequence
    report(R_fatal,"Unknown option: %s\n",a1);
    config_error++; return 1;
}

static int handle_options(int argc, const char *argv[])
{int c; int val;
    config_error=0; // no error yet
    memset(&GlobalParams,0,sizeof(GlobalParams)); 
    // some easy cases
    if(argc<2){ // no options
        short_help(); return 1;
    }
    if(strcmp(argv[1],"--version")==0){version(); return 1;}
    if(strcmp(argv[1],"--dump")==0){ dump_config(); return 1;}
    if(strcmp(argv[1],"--help=vlp")==0){ vlp_help(); return 1; }
    if(strcmp(argv[1],"--help=out")==0){ out_help(); return 1; }
    if(strcmp(argv[1],"--help")==0){ long_help(); return 1; }
    if(strcmp(argv[1],"-h")==0 || strcmp(argv[1],"-help")==0){ 
        short_help(); return 1; }
    for(c=1;c<argc;c++)if(argv[c][0]!='-'){ // vlp file
        if(PARAMS(VlpFile)){
            report(R_fatal,"Only one input file can be specified\n");
            config_error++; return -1;
        }
        PARAMS(VlpFile)=argv[c];
    } else if(argv[c][1]=='-'){ // option starts with --
        if(strncmp(argv[c],"--config=",9)==0){
            PARAMS(ConfigFile)=argv[c]+9; // check if not \0
        } else if(strncmp(argv[c],"--name=",7)==0){
            PARAMS(ProblemName)=argv[c]+7;
        } else if(strncmp(argv[c],"--boot=",7)==0){
            PARAMS(BootFile)=argv[c]+7;
        } else { // --KEYWORD=value
            int r=treat_keyword(argv[c]+2);
            if(r==-1){
                report(R_fatal,"Keyword value is out of range: %s\n",argv[c]);
                config_error++; return -1;
            }
            if(r==0){
                report(R_fatal,"Unknown keyword or invalid value: %s\n",argv[c]);
                config_error++; return -1;
            }
        }
    } else switch(argv[c][1]){ // option starts with a single '-'
        case 'c': // config file
            if(argv[c][2]){
                report(R_fatal,"Unknown option: %s\n",argv[c]);
                config_error++; return -1;
            }
            if(argv[c+1]==NULL){
                report(R_fatal,"Missing argument after option %s\n",argv[c]);
                config_error++; return -1;
            }
            c++; PARAMS(ConfigFile)=argv[c];
            break;
        case 'o': // output file
            val=argv[c][2];
            if(val && ((val!='v' && val!='f') || argv[c][3])){
                report(R_fatal,"Unknown option: %s\n",argv[c]);
                config_error++; return -1;
            }
            if(argv[c+1]==NULL){
                report(R_fatal,"Missing argument after option %s\n",argv[c]);
                config_error++; return -1;
            }
            c++; if(val=='v') PARAMS(SaveVertexFile)=argv[c];
            else if(val=='f') PARAMS(SaveFacetFile)=argv[c];
            else PARAMS(SaveFile)=argv[c];
            break;
        case 'n': // problem name
            if(argv[c][2]){
                report(R_fatal,"Unknown option: %s\n",argv[c]);
                config_error++; return -1;
            }
            if(argv[c+1]==NULL){
                report(R_fatal,"Missing argument after option: %s\n",argv[c]);
                config_error++; return -1;
            }
            c++; PARAMS(ProblemName)=argv[c];
            break;
        case 'm': // message level
            if(integer_arg(argv[c],argv[c+1],&c,&val)) return -1;
            if(val<0 || val > 3){
                report(R_fatal,"Argument of -m is out of range\n");
                config_error++; return -1;
            }
            PARAMS(ARGm)=val; PARAMS(ARGm_set)=1;
            break;
        case 'q': // quiet, same as -m0
            if(argv[c][2]){
                report(R_fatal,"Unknown option: %s\n",argv[c]);
                config_error++; return -1;
            }
            PARAMS(ARGm)=0; PARAMS(ARGm_set)=1;
            break;
        case 'p': // progress report
            if(integer_arg(argv[c],argv[c+1],&c,&val)) return -1;
            if(val<0 || val>1000000 || (0<val && val<5)){
                report(R_fatal,"Argument of -p is out of range\n");
                config_error++; return -1;
            }
            PARAMS(ARGp)=val; PARAMS(ARGp_set)=1;
            break;
        case 'y': // report vertices as generated
            val=-1;
            if(strcmp(argv[c],"-y")==0 || strcmp(argv[c],"-y+")==0
               || strcmp(argv[c],"-y1")==0) val=1;
            else if(strcmp(argv[c],"-y-")==0 
               || strcmp(argv[c],"-y0")==0) val=0;
            else {
                report(R_fatal,"Unknown option: %s\n",argv[c]);
                config_error++; return -1;
            }
            PARAMS(ARGy)= val; PARAMS(ARGy_set)=1;
            break;
        default:
            report(R_fatal,"Unknown option: %s\n",argv[c]);
            config_error++; return -1;
    }
    if(!PARAMS(VlpFile) || !*PARAMS(VlpFile)){
        report(R_fatal,"Missing input vlp file\n");
        config_error++; return -1;
    }
    read_config_file(); // parse the config file
    if(config_error) return -1;
    set_default_values(); // set default values for undefined ones
    return 0;
}
/***********************************************************************
* Process command line options and config parameters
*
* int postprocess_parameters(void)
*    make final adjustment to parameters; 
*    figure out problem name, output filenames;
*    issue warning if no output is requested, or out of memory
*
* int process_parameters(int argc, char *argv[])
*    call handle_options(), then postprocess_parameters()
*    Return the verdict: continue (0), done (1), error (-1)
*/

static void postprocess_parameters(void)
{
    if(PARAMS(ARGm_set)){ // -m[0..3]
        PARAMS(OracleMessage)=PARAMS(ARGm);
        PARAMS(MessageLevel)=PARAMS(ARGm);
        if(PARAMS(ARGm)==0) PARAMS(PrintStatistics)=0;
    }
    if(PARAMS(ARGy_set)){ // -y+ or -y-
        PARAMS(VertexReport)=PARAMS(ARGy);
    }
    if(PARAMS(ARGp_set)){ // -p T
        PARAMS(ProgressReport) = PARAMS(ARGp);
    }
    // empty strings are replaced by NULL
    if(PARAMS(BootFile) && !*PARAMS(BootFile)) PARAMS(BootFile)=0;
    if(PARAMS(SaveFile) && !*PARAMS(SaveFile)) PARAMS(SaveFile)=0;
    if(PARAMS(SaveVertexFile) && !*PARAMS(SaveVertexFile)) PARAMS(SaveVertexFile)=0;
    if(PARAMS(SaveFacetFile) && !*PARAMS(SaveFacetFile)) PARAMS(SaveFacetFile)=0;
    if(PARAMS(SaveFile)){ // -o <file>
        if(PARAMS(SaveVertexFile) && 
           strcmp(PARAMS(SaveFile),PARAMS(SaveVertexFile))==0){
            PARAMS(SaveVertexFile)=NULL;
            if(!PARAMS(SaveVertices)) PARAMS(SaveVertices)=1;
        }
        if(PARAMS(SaveFacetFile) && 
           strcmp(PARAMS(SaveFile),PARAMS(SaveFacetFile))==0){
            PARAMS(SaveFacetFile)=NULL;
            if(!PARAMS(SaveFacets)) PARAMS(SaveFacets)=1;
        }
    } else { // no -o, we might have -ov or -of 
        if(!PARAMS(SaveVertexFile)) PARAMS(SaveVertices)=0;
        if(!PARAMS(SaveFacetFile)) PARAMS(SaveFacets)=0;
    }
    // do we have any output?
    if(!PARAMS(VertexReport) && !PARAMS(PrintVertices)
       && !PARAMS(SaveVertices) && !PARAMS(SaveVertexFile)
       && !PARAMS(SaveFacets) && !PARAMS(SaveFacetFile)){
        report(R_fatal,"No output is specified; all computation would be lost...\n");
        config_error++; return;
    }
    if(PARAMS(SaveFile) && !PARAMS(SaveVertices) && !PARAMS(SaveFacets)){
        report(R_fatal,"No content is specified for the output file, it would be empty...\n");
        config_error++; return;
    }
    // create problem name
    if(PARAMS(ProblemName)==NULL){
        int l,i; const char *n=PARAMS(VlpFile);
        for(l=0,i=0;n[i];i++){
            if((n[i]=='/' || n[i]=='\\') && n[i+1]) l=i+1;
        }
        PARAMS(ProblemName)=n+l;
    }
}

int process_parameters(int argc, const char *argv[])
{int v;
    v=handle_options(argc,argv);
    if(v==0){
        postprocess_parameters();
        if(config_error) v=-1;
    }
    if(v==-1){
         report(R_fatal,"Use '" PROGNAME " --help' for a complete list of options.\n");
    }
    return v;
}

/***********************************************************************
* Show parameter values which differ from their default values
*
* void sho_parameters(void)
*    report algorithm and oracle parameters only
*/
void show_parameters(char *hdr)
{   /* integer and bool parameters */
#define CFG(x)	\
    if(PARAMS(x) != DEF_##x){ report(R_txt,"%s " #x " = %d\n",hdr,PARAMS(x)); hdr=""; }
    CFG(OracleMethod);		/* primal / dual */
    CFG(OraclePricing);
    CFG(OracleRatioTest);
    CFG(OracleItLimit);
    CFG(OracleTimeLimit);
    CFG(OracleScale);		/* scale the constraint matrix */
    CFG(ShuffleMatrix);		/* random shuffle of the constraint matrix */
    CFG(RoundVertices);		/* round vertices reported by the oracle */
    CFG(RandomFacet);		/* pick next facet randomly */
    CFG(ExactFacetEq);		/* recompute facet equation immediately */
    CFG(VertexPoolSize);	/* use vertex pool */
    CFG(OracleCallLimit);	/* oracle call limit per iteration */
    CFG(RecalculateFacets);	/* how ofter recalculate facets */
#undef CFG
    /* double parameters */
#define CFG(x)	\
    if(PARAMS(x) != DEF_##x){ report(R_txt,"%s " #x " = %lg\n",hdr,PARAMS(x)); hdr=""; }
    CFG(PolytopeEps);		/* facet and vertex adjacency */
    CFG(ScaleEps);		/* rounding when retrieving a facet equation */
    CFG(LineqEps);		/* tolerance in system of linear equations */
    CFG(RoundEps);		/* tolerance for RoundVertices */
#undef CFG
}

/* EOF */
