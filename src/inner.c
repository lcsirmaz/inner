/** inner.c inner approximation using glpk oracle **/

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "report.h"
#include "inner.h"
#include "data.h"
#include "params.h"
#include "poly.h"
#include "glp_oracle.h"

/*********************************************************************
* Initialize random; get elapsed time in 0.01 seconds, nice printing
*
* void initialize_random()
*    initialize random by asking the time and the pid of the process
*
* unsigned long timenow
*    elapsed time in 0.01 seconds, set by gettime100()
*
* unsigned long gettime100()
*    set timenow and return the elapsed time
*
* char *showtime(unsigned long t)
*    readable printout of the time in days, hours, minutes and seconds
*
* char *readable(double w, int slot)
*    print w using k,M,G,P postfix using one of the slots 0,1,2,3
*/

#include <sys/time.h> /* gettimeofday() */
#include <unistd.h>   /* pid() */
/** initialize random **/
inline static void initialize_random(void)
{struct timeval tv;
    gettimeofday(&tv,NULL);
    srandom( PARAMS(TrueRandom) ? getpid()^tv.tv_sec^tv.tv_usec : 0x12345678 );
}
/** the elapsed time in 0.01 seconds */
unsigned long timenow=0;
/** set timenow and return the elapsed time */
static unsigned long gettime100(void)
{struct timeval tv; static unsigned long starttime=0;
    if(gettimeofday(&tv,NULL)) return timenow; // some problem
    if(starttime){
        timenow = (tv.tv_sec*100 + (tv.tv_usec+5000u)/10000u)-starttime;
        return timenow;
    }
    starttime=tv.tv_sec*100 + (tv.tv_usec+5000u)/10000u;
    timenow = 0;
    return timenow;
}
/** print time in pretty format **/
static char *showtime(unsigned long t)
{static char buff[50]; int m,s;
    if(t<60*100){ sprintf(buff,"%.2f",((double)t)*0.01); return buff;}
    s=(int)((t+50)/100); m=s/60; s=s%60;
    if(m<60){ sprintf(buff,"%d:%02d",m,s); return buff; }
    if(m<60*24){ sprintf(buff,"%d:%02d:%02d",m/60,m%60,s); return buff;}
    sprintf(buff,"%dd%02d:%02d:%02d",m/1440,(m/60)%24,m%60,s); return buff;
}
/** print a huge double with postfix k,M,G,P to a slot 0..3 **/
static char* readable(double w, int slot)
{static char slots[4][30]; char *buff;
    buff=slots[slot];
    if(w<0.0){ w=0.0; } /* only >=0 numbers */
    if(w<1000.0){ sprintf(buff,"%.2lf",w); return buff;}
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfk",w); return buff; }
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfM",w); return buff; }
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfG",w); return buff; }
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfP",w); return buff; }
    /* too large value */
    sprintf(buff,"%lgP",w); return buff;
}

/*******************************************************************
* Timing the oracle calls
*
* unsigned long oracletime
*    total time required by the oracle in 0.01 seconds
* int oraclecalls
*    number of times the oracle was called
* int ask_oracle_with_timer()
*    ask oracle, measure the elapsed time, print error messages
*/
static unsigned long oracletime=0;
static int oraclecalls=0;

inline static int ask_oracle_with_timer(void)
{unsigned long thistime; int res;
    oraclecalls++;
    thistime=gettime100();
    res=ask_oracle();
    oracletime += gettime100()-thistime;
    return res;
}

/*******************************************************************
* Progress report
*
* unsigned long progresstime, progressdelay
*    when the last progress report was made
*    minimum delay between two reports, in 0.01 seconds
*
* unsignel long chktime, chkdelay
*    when the last checkpoint was created, and the minimum delay
*    between two checkpoints
*
* int poolstat
*    number of entries in VertexPool, set by find_next_vertex()
*
* int facetstat
*    whether new facet stat is available after calling add_new_vertex()
*
* void progress_stat(void)
*    show progress report, save the report time
*
* void progress_stat_if_expired(int ask_time)
*    show progress report when expired
*
* void report_new_vertex(int ask_time)
*    report coordinates of the new vertex when requested; give 
*    status report when the delay has passed
*
* void report_new_facet(int facetno)
*    report facet equation for the new final facet; give status
*    report when the delay has passed
*
* void report_memory(void)
*    when total memory allocated changes, call report_memory_usage()
*
* int memory_limit_reached(void)
*    return 1 if the allocated memory exceeds the one set in PARAMS(MemoryLimit)
*/
static unsigned long progresstime=0, progressdelay=0;
static unsigned long chktime=0, chkdelay=0;

static int poolstat=0, facetstat=0;

static void progress_stat(void)
{   progresstime=timenow;
    get_dd_facetno();
    report(R_txt,"I%8.2f] Elapsed: %s, vertices: %d, facets final: %d, pending: %d",
        0.01*(double)timenow,
        showtime(timenow),
        dd_stats.vertexno,
        dd_stats.final_facets_no,
        dd_stats.living_facets_no-dd_stats.final_facets_no);
    if(poolstat){ report(R_txt,", pool: %d",poolstat); poolstat=0; }
    if(facetstat){
        report(R_info,", eq: %d, out: %d, in: %d",
           dd_stats.facet_zero, dd_stats.facet_neg,dd_stats.facet_new);
        facetstat=0;
    }
    report(R_txt,"\n");
}

static void progress_stat_if_expired(int ask_time)
{   if(PARAMS(ProgressReport)==0) return;
    if(ask_time) gettime100();
    if(timenow-progresstime < progressdelay) return;
    progress_stat();
    flush_report();
}

static void report_new_vertex(int ask_time)
{int flush=0;
    if(PARAMS(ProgressReport)==0 && !PARAMS(VertexReport)) return;
    if(PARAMS(ProgressReport)){
      if(ask_time){ gettime100(); ask_time=0; }
      if(timenow-progresstime >= progressdelay){
        progress_stat();
        flush=1;
      }
    }
    if(PARAMS(VertexReport)){
      if(ask_time) gettime100();
      report(R_txt,"[%8.2f] ",0.01*(double)timenow);
      print_vertex(R_txt,(PARAMS(Direction) ? -1.0 : +1.0),VertexOracleData.overtex);
      flush=1;
    }
    if(flush) flush_report();
}

static void report_new_facet(int fno)
{int flush=0;
    if(PARAMS(ProgressReport)==0 && PARAMS(FacetReport)==0) return;
    gettime100();
    if(PARAMS(ProgressReport)){
       if(timenow-progresstime >= progressdelay){
         progress_stat(); flush=1;
       }
    }
    if(PARAMS(FacetReport)){
      report(R_txt,"[%8.2f] ",0.01*(double)timenow);
      print_facet(R_txt,fno); flush=1;
    }
    if(flush) flush_report();
}

inline static void report_memory(void)
{static int last_memreport=0;
    if(! PARAMS(MemoryReport) || 
       last_memreport==dd_stats.memory_allocated_no ) return;
    last_memreport=dd_stats.memory_allocated_no;
    report(R_txt,"M%8.2f] ", 0.01*(double)timenow);
    report_memory_usage();
    flush_report();
}

inline static int memory_limit_reached(void)
{  return PARAMS(MemoryLimit)>=100 && 
      (((double)dd_stats.total_memory)*0.000001>((double)PARAMS(MemoryLimit)))
   ? 1 : 0;
}

/*******************************************************************
* Dump and save the result, print statistics
*
* EQSEP, DASHSEP
*    separators made of = and - respectively
*
* void dump_and_save(int how)
*    dump and save final vertices and facets and print statistics.
*    The argument tells how the program terminated:
*      0 - normal
*      1 - error but data is consistent
*      2 - other error (numerical or out of memory)
*      3 - interrupt. 
*/

static void print_vertexpool_content(report_type rt); /* forward declaration */

#define EQSEP   "================================"
#define DASHSEP "--------------------------------"

static void dump_and_save(int how)
{unsigned long endtime; int partial;
    endtime=gettime100(); // program finished
    if(PARAMS(ProgressReport)) progress_stat();
    partial = how==0 ? 0 : 1; // print data when completed
    if(PARAMS(PrintVertices) > partial){
        if(partial) report(R_txt,"Partial list of vertices:\n");
        report(R_txt, "\n" DASHSEP "\n");
        print_vertices(R_txt);
        if(partial && PARAMS(VertexPoolSize)>=5 ){
            report(R_txt,"\nAdditional vertices in the pool:\n");
            print_vertexpool_content(R_txt);
        }
    }
    if(PARAMS(PrintFacets) > partial){
        if(partial) report(R_txt,"Partial list of facets:\n");
        report(R_txt, "\n" DASHSEP "\n");
        print_facets(R_txt);
    }
    if(PARAMS(PrintStatistics)){ /* statistics, only if not quiet */
      report(R_txt,"\n" DASHSEP "\n"
      "Problem %s\n"
      " name                    %s\n"
      " output                  %s\n"
      "%s%s%s"
      "%s%s%s"
      " rows, cols, objs        %d, %d, %d\n"
      " vertices, facets        %d, %d\n"
      " total time              %s\n",
      how==0 ? "completed" : how<=2 ? "aborted with error" : "interrupted",
      PARAMS(ProblemName),
      PARAMS(SaveFile) ? PARAMS(SaveFile) : 
        PARAMS(SaveVertexFile)||PARAMS(SaveFacetFile) ? "" : "[none]",
      PARAMS(SaveVertexFile) ? "   vertices              " : "",
      PARAMS(SaveVertexFile) ? PARAMS(SaveVertexFile) : "",
      PARAMS(SaveVertexFile) ? "\n": "",
      PARAMS(SaveFacetFile) ?  "   facets                " : "", 
      PARAMS(SaveFacetFile) ? PARAMS(SaveFacetFile) : "",
      PARAMS(SaveFacetFile) ? "\n" : "",
      PARAMS(ProblemRows), PARAMS(ProblemColumns), PARAMS(ProblemObjects),
      vertex_num(), facet_num(), showtime(endtime));
      report(R_txt, DASHSEP "\nStatistics\n"
      "LP:\n"
      " oracle calls            %d\n"
      "   avg iterations/calls  %s\n"
      " total oracle time       %s\n",
      oraclecalls, readable((0.0001+get_oracle_rounds())/(0.0001+oraclecalls),0),
      showtime(oracletime));
      report(R_txt,
      "Combinatorics:\n"
      " vertices added          %d\n"
      " facets probed           %d\n"
      " facets # max            %d\n"
      " memory allocated        %s%s\n",
      dd_stats.iterations+1,dd_stats.facetenquires,
      dd_stats.max_facets,
      readable(dd_stats.total_memory,0),
      dd_stats.out_of_memory ? " (out of memory)" : "");
      if(dd_stats.instability_warning) report(R_txt,
      " instability warnings    %d\n",
      dd_stats.instability_warning);
      report(R_txt,
      " storage expansion\n"
      "   vertices # / upto     %d / %d\n"
      "   facets   # / upto     %d / %d\n"
      " facet manipulating\n"
      "   added avg / max       %s / %s\n"
      "   compressed            %d\n"
      " number of ridge tests\n"
      "   avg / max             %s / %s\n",
      dd_stats.vertices_allocated_no,dd_stats.vertices_allocated,
      dd_stats.facets_allocated_no,dd_stats.facets_allocated,
      readable(dd_stats.avg_facetsadded,0),readable(dd_stats.max_facetsadded,1),
      dd_stats.facets_compressed_no,
      readable(dd_stats.avg_tests,2),readable(dd_stats.max_tests,3));
    }
    if(PARAMS(PrintParams)){
        show_parameters(DASHSEP "\nParameters with non-default values\n");
    }
    partial = how<2 ? 0 : 1; // save data when consistent
    if(PARAMS(SaveVertices)>partial || (partial==0 && PARAMS(SaveVertexFile))){
        report(R_savevertex,"C name=%s, rows=%d, cols=%d, objs=%d\n"
          "C vertices=%d, facets=%d\n\n",
          PARAMS(ProblemName), 
          PARAMS(ProblemRows), PARAMS(ProblemColumns), PARAMS(ProblemObjects),
          vertex_num(), facet_num());
        if(how) report(R_savevertex,"C *** Partial list of vertices ***\n");
        print_vertices(R_savevertex);
        if(how && PARAMS(VertexPoolSize)>=5 ){
            report(R_savevertex,"\nC *** Additional vertices in the pool ***\n");
            print_vertexpool_content(R_savevertex);
        }
        report(R_savevertex,"\n");
    }
    if(PARAMS(SaveFacets)>partial || (partial==0 && PARAMS(SaveFacetFile))){
        report(R_savefacet,"C name=%s, rows=%d, cols=%d, objs=%d\n"
          "C vertices=%d, facets=%d\n\n",
          PARAMS(ProblemName), 
          PARAMS(ProblemRows), PARAMS(ProblemColumns), PARAMS(ProblemObjects),
          vertex_num(), facet_num());
        if(how) report(R_savefacet,"C *** Partial list of facets ***\n");
        print_facets(R_savefacet);
        report(R_savefacet,"\n");
    }
    close_savefiles();
}

/**************************************************************************
* Find the next vertex to be added
*
* int VertexPoolAfter = 20
* int VertexPoolMinFacets = 1000
*   use vertexpool only when we have at least that many vertices. or
*   that many unprocessed facets
*   
* vertexpool_t vertexpool[VertexPoolSize]
*   vertices known but not added yet to the approximation.
*
* int facets_recalculated
*   TRUE just after calling recalculate_facets(); and set to FALSE
*   when calling add_new_vertex()
*
* int init_vertexpool()
*   allocates memory to the vertex pool. Returns non-zero if out of
*   memory.
*
* void print_vertexpool_content(void)
*   report vertices in the vertex pool.
*
* int same_vec(double v1[0:DIM-1], double v2[0:DIM-1])
*   checks whether the two vectors are the same (1) or not (0)
*
* int next_vertex_coords(int checkVertexPool)
*   try to generate a new vertex. Ask the oracle about the next
*   facet returned by get_next_facet(-1). Check whether facet has
*   been asked before if checkVertexPool is set. If the returned
*   vertex is on the facet, mark the facet as final and repeat.
*   Return value:
*     0:  there are no more facets
*     1:  algorithm interrupted
*     2:  unbounded
*     4:  some error (oracle failed, computational error, etc)
*     5:  facet was encountered before
*     6:  next vertex is VertexOracleData.overtex[]
*     7:  memory limit is hit
*
* int fill_vertexpool(int limit)
*   fill the vertex pool, limiting unsuccessful oracle calles to limit.
*   Return value:
*     0:  pool is filled (maybe empty if no more vertex)
*     1:  interrupt
*     2:  problem unbounded
*     4:  numerical error
*     7:  memory limit is hit
*
* int find_next_vertex(void)
*   finds the next vertex to be added to the approximating polytope.
*   Without the vertex pool it returns next_vertex_coords().
*   Otherwise returns the vertex with the largest probe_vertex(v)
*   value. Returned values are the same as for next_vertex_coords().
*/

#define VertexPoolAfter	20  /* use vertexpool only after that many vertices */
#define VertexPoolMinFacets 1000 /* or after that many unprocessed facets */

typedef struct {
    int    occupied;	/* 0=no, 1=yes */
    double *coords;	/* pointer to vertex coordinates */
    double *facet;	/* pointer to the facet coordinates */
} vertexpool_t;

static vertexpool_t *vertexpool=NULL;

#define DIM	PARAMS(ProblemObjects) /* problem dimension */

static int facets_recalculated=0;

static int init_vertexpool(void) /* call only when DIM has been set */
{int i; double *pool;
    if(PARAMS(VertexPoolSize)<5) return 0; // don't use it
    vertexpool=malloc(PARAMS(VertexPoolSize)*sizeof(vertexpool_t));
    pool=malloc(PARAMS(VertexPoolSize)*2*DIM*sizeof(double));
    if(!vertexpool || !pool){
        report(R_fatal,"init_vertexpool: out of memory\n");
        return 1;
    }
    for(i=0;i<PARAMS(VertexPoolSize);i++){
        vertexpool[i].occupied=0;
        vertexpool[i].coords=pool; pool+=DIM;
        vertexpool[i].facet=pool;  pool+=DIM;
    }
    return 0;
}

static void print_vertexpool_content(report_type where)
{int i; double dir;
    dir= PARAMS(Direction) ? -1.0 : 1.0;
    for(i=0;i<PARAMS(VertexPoolSize);i++) if(vertexpool[i].occupied){
        report(where,"v ");
        print_vertex(where,dir,vertexpool[i].coords);
    }
}

inline static int same_vec(double f1[], double f2[])
{int i; double d;
    for(i=0;i<DIM;i++){
       d=f1[i]-f2[i];
        if(d> PARAMS(PolytopeEps) || d < -PARAMS(PolytopeEps) )
            return 0; /* no */
    }
    return 1; /* yes */
}

static int next_vertex_coords(int checkVertexPool)
{int i,j; double d,dd; int boottype; int pos_facet=0;
again:
    if(dobreak) return 1; /* break meanwhile */
    if(memory_limit_reached()) return 7; /* memory limit reached */
    // boot file
    while(nextline(&boottype)) if(boottype==1){ // V line
        if(parseline(PARAMS(ProblemObjects),VertexOracleData.overtex)){
            return 4; // error
        }
        memset(VertexOracleData.ofacet,0,DIM*sizeof(double));
        return 6;
    }
    j=get_next_facet(-1);
    if(j<0) return 0; /* terminated successfully */
    get_facet_into(j,VertexOracleData.ofacet);
    if(checkVertexPool){ /* ask oracle only when not asked before */
        for(i=0;i<PARAMS(VertexPoolSize);i++) if(vertexpool[i].occupied
           && same_vec(vertexpool[i].facet,VertexOracleData.ofacet)){
            return 5;
        }
    }
    if((i=ask_oracle_with_timer())!=0){ // error
        return i==ORACLE_UNBND ? 2 : 4;
    }
    d=VertexOracleData.ofacet[DIM]; dd=0.0;
    for(i=0;i<DIM;i++){
        dd+= VertexOracleData.ofacet[i];
        d+= VertexOracleData.ofacet[i]*VertexOracleData.overtex[i];
    }
    d /= dd; // sum of facet coeffs is 1.0 
    if(PARAMS(PolytopeEps) < d){ /* numerical error */
        // make sure facets are recalculated immediately before issuing an error
        if(facets_recalculated){
            pos_facet++;
            if(pos_facet<3){
                report(R_warn,"New vertex is on positive side of facet %d (d=%lg), trying another facet ...\n",j,d);
                dd_stats.instability_warning++;
                goto again;
            }
            report(R_fatal,"Numerical error: new vertex is on the positive side of facet %d (d=%lg)\n",j,d);
            return 4;
        }
        report(R_warn,"New vertex is on positive side of facet %d (d=%lg), recalculating facets ...\n",j,d);
        dd_stats.instability_warning++;
        recalculate_facets();
        if(dd_stats.out_of_memory || dd_stats.numerical_error){
            return 4; // error during computation
        }
        facets_recalculated=1;
        goto again;
    }
    if(-PARAMS(PolytopeEps) < d){ /* vertex is on the facet */
        mark_facet_as_final(j);
        report_new_facet(j);// progress_stat_if_expired(1);
        goto again;
    }
    return 6;
}

static int fill_vertexpool(int limit)
{int i,ii; int oracle_calls=0;
    for(i=0;i<PARAMS(VertexPoolSize);i++)if(!vertexpool[i].occupied){
        switch(next_vertex_coords(1)){
      case 0: return 0; /* no more vertices, done */
      case 1: return 1; /* break */
      case 2: return 2; /* unbounded */
      case 4: return 4; /* error */
      case 7: return 7; /* memory limit */
      case 5: break;    /* same facet encountered again, skip it */
      default:
              for(ii=0;ii<PARAMS(VertexPoolSize);ii++) if(vertexpool[ii].occupied
                 && same_vec(VertexOracleData.overtex,vertexpool[ii].coords)) break;
              if(ii==PARAMS(VertexPoolSize)){
                  memcpy(vertexpool[i].coords,VertexOracleData.overtex,DIM*sizeof(double));
                  memcpy(vertexpool[i].facet,VertexOracleData.ofacet,DIM*sizeof(double));
                  vertexpool[i].occupied=1;
              } else { // got the same vertex; save the newer values
                  memcpy(vertexpool[ii].coords,VertexOracleData.overtex,DIM*sizeof(double));
                  memcpy(vertexpool[ii].facet,VertexOracleData.ofacet,DIM*sizeof(double));
                  // and check how many unsuccessful calls we have made
                  oracle_calls++;
                  if(limit && oracle_calls>=limit) return 0; // done
              }
    }}
    return 0;
}

static int find_next_vertex(void)
{int i,maxi,cnt; int w,maxw;
    if(PARAMS(VertexPoolSize)<5 || ( dd_stats.vertexno < VertexPoolAfter &&
         dd_stats.facet_zero+dd_stats.facet_pos+dd_stats.facet_new < VertexPoolMinFacets)
       )
        return next_vertex_coords(0);
    // fill the vertex pool
    i=fill_vertexpool(PARAMS(OracleCallLimit));
    if(i) return i; // break or numerical error
    /* find the weight of stored vertices */
    maxi=-1; maxw=0; cnt=0;
    for(i=0;i<PARAMS(VertexPoolSize);i++)if(vertexpool[i].occupied){
        cnt++;
        w=probe_vertex(vertexpool[i].coords);
        if(maxi<0 || maxw<w){ maxi=i; maxw=w; }
    }
    poolstat=cnt;
    if(maxi<0) return 0; // no more vertices
    vertexpool[maxi].occupied=0;
    memcpy(VertexOracleData.overtex,vertexpool[maxi].coords,DIM*sizeof(double));
    return 6; // next vertex is in VertexOracleData.overtex
}

#undef DIM
/**************************************************************************
* The main loop of the algorithm
* int inner(void)
*   when it starts, all parameters in PARAMS have been set. The steps:
*   o  read_vlp() reads in the the MOLP problem form a vlp file
*   o  check that output files are writable
*   o  set_oracle_parameters() sets the oracle parameters
*   o  the first direction to probe has all coordinates set to 1.0
*   o   ask_oracle() returns the first vertex of the final polytope
*   o   init_dd() creates the first approximating simplex using this
*         vertex and the positive ideal endpoints of the coordinate axes
*   loop:
*   o  get_next_facet(-1) returns the next facet of the approximating
*        polytope which is not known to be final. If none, the algorithm
*        terminates. Otherwise
*   o  ask_oracle() returns an extremal vertex at the direction of the
*        facet. If the vertex is on the facet, then
*   o  mark_facet_as_final(), and goto loop.
*        If the vertex is on the negative side of the facet then
*   o  add_new_vertex(), and goto loop.
*
* int break_inner(int how)
*   when the inner routine is interrupted by the signal, this procedure
*   kicks in. It goes over all facets of the actual approximation,
*   and calls the oracle to produce potentially new vertices. Argument
*   tells if interrupt or memory limit reached.
*
* int handle_new_vertex(void)
*   the new vertex is in VertexOracleData.overtex; make reports, add as
*   a new vertex by calling add_new_vertex(), take care of timed actions
*   such as reports and checkpoints.
*   Return values:
*     1:  OK
*     0:  error during computation (cannot continue)
*/

/** the main loop was interrupted; return 3, 4, 5 **/
static int break_inner(int how)
{int j; unsigned long aborttime;
    aborttime=gettime100(); dobreak=0;
    report(R_fatal,"\n\n" EQSEP "\n %s after %s, vertices=%d, facets=%d\n",
      how==0 ? "Program run was interrupted" : "Memory limit reached",
      showtime(aborttime), vertex_num(), facet_num());
    if(!PARAMS(ExtractAfterBreak)) return 5; // normal termination
    // don't do if no need to extract data
    if(PARAMS(PrintVertices)<2 && PARAMS(SaveVertices)<2 &&
       !PARAMS(VertexReport)){
        report(R_fatal,"Result of postprocessing would be lost, not doing...\n");
        return 5;
    }
    report(R_fatal,"Checking additional vertices. This may take some time...\n"
      EQSEP "\n\n");
    // at least one of PrintVertices and SaveVertices should be set
    if(PARAMS(PrintVertices)<2 && PARAMS(SaveVertices)<2)
        PARAMS(PrintVertices)=2;
    /* get vertices from the pool */
    if(PARAMS(VertexPoolSize)>=5){
        for(j=0;j<PARAMS(VertexPoolSize);j++) if(vertexpool[j].occupied
            && store_vertex(vertexpool[j].coords)){
            memcpy(VertexOracleData.overtex,vertexpool[j].coords,PARAMS(ProblemObjects)*sizeof(double));
            report_new_vertex(1);
        }
    }
    /* search for new vertices by calling the oracle for all facets */
    j=-1;
    while((j=get_next_facet(j+1))>=0){
        if(dobreak){
            dobreak=0;
            if(gettime100()-aborttime>50){
                report(R_fatal,"\n" EQSEP "\n"
                  "Post-processing was interrupted after %s\n",
                  showtime(timenow-aborttime));
                return 7; // postprocess aborted
            }
        }
        get_facet_into(j,VertexOracleData.ofacet);
        mark_facet_as_final(j);
        if(ask_oracle_with_timer()){ // error
            return 6; // error during postprocess
        }
        if(store_vertex(VertexOracleData.overtex))
            report_new_vertex(0);
        else
            progress_stat_if_expired(0);
        if(dd_stats.out_of_memory && !PARAMS(VertexReport)){
            return 6; // error in postprocess
        }
    }
    return dd_stats.out_of_memory ? 6 : 5; // terminated normally
}

static int handle_new_vertex(void)
{   // progress report plus info about the new vertex
    report_new_vertex(0);
    add_new_vertex(VertexOracleData.overtex);
    gettime100(); 
    facetstat=1; progress_stat_if_expired(0);
    if(dd_stats.out_of_memory || dd_stats.numerical_error){
        return 0; // error during computation
    }
    facets_recalculated=0;
    // recalculate facets if instructed so
    if(PARAMS(RecalculateFacets)>=5 &&
       ((1+dd_stats.iterations)%PARAMS(RecalculateFacets))==0){
        // report what we are going to do
        report(R_info,"I%8.2f] recalculating facets...\n",
             0.01*(double)timenow);
        recalculate_facets();
        gettime100();
        if(dd_stats.out_of_memory || dd_stats.numerical_error){
            return 0; // error during computation
        }
        facets_recalculated=1;
    }
    // check consistency if instructed
    if(PARAMS(CheckConsistency)>=5 &&
       ((1+dd_stats.iterations)%PARAMS(CheckConsistency))==0){
        // report what we are going to do
        report(R_warn,"I%8.2f] checking data consistency...\n",
              0.01*(double)timenow);
        if(check_consistency()) { // error
            report(R_fatal,"Consistency error: data structure has numerical errors.\n");
            return 0;
        }
        gettime100();
    }
    report_memory();
    return 1;
}

/** the main algorithm **/
typedef enum {
inp_none,	/* no special input is given */
inp_boot,	/* reading verives from --boot */
inp_resume	/* vertices and facets form --resume */
} input_type_t;

int inner(void)
{int i; int retvalue=0; input_type_t inp_type=inp_none;
    initialize_random(); // initialize random numbers
    if(read_vlp()){ return 1; } // data error before start
    if(check_outfiles()){ return 1;}
    if(PARAMS(BootFile)){ // we have a bootfile
        if(init_reading(PARAMS(BootFile))){ return 1; } // error
        inp_type=inp_boot;
    } else if(PARAMS(ResumeFile)){ // we have a resume file
        if(init_reading(PARAMS(ResumeFile))){ return 1; }
        inp_type=inp_resume;
    }
    report(R_info,"C MOLP problem=%s, %s\n"
       "C rows=%d, columns=%d, objectives=%d\n",
       PARAMS(ProblemName),PARAMS(Direction)?"maximize":"minimize",
       PARAMS(ProblemRows),PARAMS(ProblemColumns),PARAMS(ProblemObjects));
    set_oracle_parameters();
    gettime100(); // initialize elapsed time, reading vlp not included
    if(init_vertexpool()){ return 1; } // fatal error
#ifdef USETHREADS
    if(create_threads()) return 1;
#endif
    progressdelay = 100*(unsigned long)PARAMS(ProgressReport);
    chkdelay = 100*(unsigned long)PARAMS(CheckPoint);
#define DIM	PARAMS(ProblemObjects) /* problem dimension */
    if(inp_type==inp_boot){ // read first vertex from bootfile
       int linetype=0;
       while(nextline(&linetype) && linetype!=1); // next V line
       if(linetype!=1){ // error
           report(R_fatal,"No vertex line was found in boot file %s\n",PARAMS(BootFile));
           retvalue=1; goto leave;
       }
       if(parseline(DIM,VertexOracleData.overtex)){
           retvalue=1; goto leave; // error
       }
       report_new_vertex(1); // take care of reporting
       if(init_dd(DIM,VertexOracleData.overtex)){ retvalue=1; goto leave; }
    } else if(inp_type==inp_resume){ // resume an earlier computation
       int linetype=0; double args[5];
       if(! nextline(&linetype) || linetype!= 4 || parseline(5,&args[0]) ){
           report(R_fatal,"Argument line is missing from the resume file %s\n",PARAMS(ResumeFile));
           retvalue=1; goto leave;
       }
       // rows, columns, objects
       if(args[2]!=(double)PARAMS(ProblemRows) ||
          args[3]!=(double)PARAMS(ProblemColumns) ||
          args[4]!=(double)PARAMS(ProblemObjects)) {
           report(R_fatal,"Resume file %s belongs to a different MOLP problem\n",PARAMS(ResumeFile));
           retvalue=1; goto leave; 
       }
       if(init_dd_structure(DIM,(int)args[0],(int)args[1])){
           retvalue=1; goto leave;
       }
       // read the resume file
       while(nextline(&linetype)) switch(linetype){
         case 1:   // V line
            if(parseline(DIM,VertexOracleData.overtex) || 
               initial_vertex(VertexOracleData.overtex)){
                retvalue=1; goto leave; // error
            }
            break;
         case 2:   // F line
            if(parseline(DIM+1,VertexOracleData.ofacet) ||
               initial_facet(1,VertexOracleData.ofacet)){
                retvalue=1; goto leave; // error
            }
            break;
         case 3:   // f line
            if(parseline(DIM+1,VertexOracleData.ofacet) ||
               initial_facet(0,VertexOracleData.ofacet)){
                retvalue=1; goto leave; // error
            }
            break;
         default:  // ignore
            break;
       }
       gettime100(); // set elapsed time
       if(PARAMS(ProgressReport)){ progress_stat(); flush_report(); }
       if(PARAMS(VertexPoolSize)>=5 && fill_vertexpool(0)){
          report(R_fatal,"Error while filling the vertex pool\n");
          retvalue=4; goto leave;
       }
    } else { // inp_type==inp_none
       // first vertex, all directions are 1.0
       for(i=0;i<DIM;i++) VertexOracleData.ofacet[i]=1.0;
       if((i=ask_oracle_with_timer())!=0){
           report(R_fatal,"Error getting the first optimal solution.\n");
           retvalue= i==ORACLE_UNBND ? 2:
                     i==ORACLE_EMPTY ? 3: 4;
           goto leave;
       }
       report_new_vertex(0); // take care of reporting
       if(init_dd(DIM,VertexOracleData.overtex)){ retvalue=4; goto leave; }
    }
    if(DIM<2){
        PARAMS(PrintFacets)=0;
        PARAMS(SaveFacets)=0;
        dump_and_save(0); // done
        retvalue=0;
        goto leave;
    }
    chktime=gettime100(); // last checktime
again:
    gettime100(); // set elapsed time
    if(dodump){ // the data is consistent, request for a dump
        dodump=0;
        report(R_info,"I%08.2f] Dumping vertices and facets\n",0.01*(double)timenow);
        make_dump();
    }
    // make checkpoint if expired
    if(PARAMS(CheckPointStub) && chktime+chkdelay<=timenow){
        chktime=timenow;
        report(R_info,"I%08.2f] Checkpoint: dumping vertices and facets\n",0.01*(double)chktime);
        make_checkpoint(); 
    }
    switch(find_next_vertex()){
      case 0: // no more vertices
        dump_and_save(0); // done
        retvalue=0; goto leave; // terminated normally
      case 1: // break
        retvalue=break_inner(0);
        dump_and_save(3);
        goto leave;
      case 2: // problem unbounded, can be numerical error as well
        dump_and_save(2);
        retvalue=2; goto leave;
      case 4: // numerical error
        dump_and_save(2);
        retvalue=4; goto leave; // error during computation
      case 7: // memory limit
        retvalue=break_inner(1);
        dump_and_save(3);
        goto leave;
      default: // next vertex returned
        if(handle_new_vertex()) goto again;
        dump_and_save(dd_stats.data_is_consistent? 1 : 2);
        retvalue=4; goto leave;
    }
#undef DIM
leave:
#ifdef USETHREADS
    stop_threads();
#endif
    return retvalue;
}

/* EOF */

