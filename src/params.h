/* params.h -- global parameters */

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

/***********************************************************************
* Program parameters
*
* struct params_t GlobalParams
*    global parameters for the problem solver. They define algorithm
*    behavior, how the result is saved, LP parameters, etc. Parameters
*    starting with ARG store command line arguments.
*/
struct params_t {
  char			/* bool parameters, or parameters with small values */
    MessageLevel,	/* 0: quiet, 1: error, 2: warning, 3: verbose */
    PrintParams,		/* print parameters which differ from their default values */
    PrintStatistics,	/* print statistics at the end */
    VertexAsFraction,	/* print vertex coordinates as fractions */
    VertexReport,	/* print out vertices when they are found */
    FacetReport,	/* print out final facets when found */
    MemoryReport,	/* print combinatorial memory usage  whenever changes */
    PrintVertices,	/* dump vertices at the end */
    PrintFacets,	/* dump facets at the end */
    SaveVertices,	/* save vertices at the end */
    SaveFacets,		/* save facets at the end */
    RandomFacet,	/* pick the facet to be tested randomly */
    ExactFacetEq,	/* always calculate the facet equation from adjacency vertices */
    ExtractAfterBreak,	/* continue after break with extracting vertices */
    ShuffleMatrix,	/* (oracle) shuffle rows, columns, and objective order.
			   helps numerical stability */
    RoundVertices,	/* (oracle) round vertex coordinates to the nearest rational */
    OracleMessage,	/* 0: quiet, 1: error, 2: on, 3: verbose */
    OracleScale,	/* scale constraint matrix;  0: no, 1: yes */
    OracleMethod,	/* 0: primal, 1: dual */
    OracleRatioTest,	/* 0: standard, 1: Harris */
    OraclePricing,	/* 0: standard, 1: steepest edge */
    Direction,		/* 0: minimize, 1: maximize, set by the Oracle */
    TrueRandom,		/* whether use true or deterministic generator */
    ARGm,		/* -m[0-3] option */
    ARGy,		/* -y+ or -y- */
    ARGm_set,ARGp_set,ARGy_set;
			/* whether these options were used */

  int			/* integer parameters */
    ProgressReport,	/* progress frequency in seconds, 0 means no report */
    RecalculateFacets,	/* after that many iterations do it */
    CheckConsistency,	/* after that many iterations do it */
    VertexPoolSize,	/* vertex pool size, 0 means no vertex pool */
    CheckPoint,		/* frequency in seconds of dumping vertices and facets */
    MemoryLimit,	/* stop when reaching that mamory usage, in Mbytes */
    TimeLimit,		/* stop when running for that many seconds */
    Threads,		/* number of threads to use, only when USETHREADS defined */
    OracleItLimit,	/* iteration limit, >=1000; =0: unlimited */
    OracleTimeLimit,	/* time limit in seconds, >=5; =0: unlimited */
    OracleCallLimit,	/* limit of oracle calls in each iteration */
    ProblemColumns,	/* problem columns, set by the Oracle */
    ProblemRows,	/* problem rows, set by the Oracle */
    ProblemObjects,	/* problem objects (dimension), set by the Oracle */
    ARGp;		/* -r N */

  double		/* double parameters */
    RoundEps,		/* (oracle) round to a rational if closer than this value, 1e-9 */
    ScaleEps,		/* coeffs in facet equation are rounded to the nearest
			   integer if they are closer than this, 3e-9 */
    PolytopeEps,	/* max distance between vertex and facet, 1.3e-8 */
    LineqEps,		/* solving linear equation for facet, 6.0*PolytopeEps */
    FacetRecalcEps;	/* report numerical instability when after recomputing
			   facet equation the old and new coords differ by
			   that much (1e-6) */

  const char		/* string parameters */
    *VlpFile,		/* the input vlp file */
    *BootFile,		/* initial list of vertices */
    *ResumeFile,	/* resume from this checkpoint file */
    *ProblemName,	/* the problem name, typically the base of the vlp file */
    *ConfigFile,	/* configuration file name */
    *CheckPointStub,	/* -oc <stub> option */
    *SaveFile,		/* -o <file> option */
    *SaveVertexFile,	/* -ov <file> option */
    *SaveFacetFile;	/* -of <file> option */
};

extern struct params_t GlobalParams;

#define PARAMS(field)	GlobalParams.field

/***********************************************************************
* Parse command line options and config file
*
* int process_parameters(int argc, char *argv[])
*    go over all command line options; handle --help, --version, --dump
*    otherwise read config file (if specified), and then
*    set default values for unset PARAMS().
*    Return value:
*      0:  options are OK, can continue
*      1:  option handled, exit normally
*     -1:  some error, error message issued, exit with error
*
* void show_parameters(char *hdr)
*    print algorithm and oracle parameters which differ from their
*    default values. Put hdr before the first line.
*/

int process_parameters(int argc, const char *argv[]);
void show_parameters(char *hdr);

/***********************************************************************
* mkstringof() macro
*    creates a string from its only argument
*/
#define stringify(x)	#x
#define mkstringof(x)	stringify(x)

/***********************************************************************
* Signals which stops processing and instruct special continuation
*/

#ifndef INNER_SIGNAL	/* interrupt computation */
#define INNER_SIGNAL	SIGUSR1
#endif
#ifndef DUMP_SIGNAL	/* dump vertices and facets */
#define DUMP_SIGNAL	SIGUSR2
#endif

/************************************************************************
* Maximal values for some parameters
*/
#ifndef MAX_THREADS
#define MAX_THREADS	64 	/* number of threads allowed */
#endif
#ifndef MAX_VERTEX_POOL
#define MAX_VERTEX_POOL	3000	/* maximum size of the vertex pool */
#endif
#ifndef MAX_OCALL_LIMIT
#define MAX_OCALL_LIMIT	100	/* unsuccessfull oracle calls per iteration */
#endif
/* EOF */

