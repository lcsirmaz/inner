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
    ReportLevel,	/* 0: quiet, 1: error, 2: warning, 3: verbose */
    PrintAsFraction,	/* print vertex coordinates as fractions */
    ShowVertices,	/* print out vertices when they are found */
    ReportMemory,	/* report combinatorial memory usage */
    DumpVertices,	/* dump vertices at the end */
    DumpFacets,		/* dump facets at the end */
    SaveVertices,	/* save vertices at the end */
    SaveFacets,		/* save facets at the end */
    RandomFacet,	/* pick the facet to be tested randomly */
    ExactFacetEq,	/* always calculate the facet equation from adjacency vertices */
    ExtractAfterBreak,	/* continue after Ctrl+C with extracting vertices */
    ShuffleMatrix,	/* (oracle) shuffle rows, columns, and objective order.
			   helps numerical stability */
    RoundVertices,	/* (oracle) round vertex coordinates to the nearest rational */
    OracleMessage,	/* 0: quiet, 1: error, 2: on, 3: verbose */
    OracleScale,	/* scale constraint matrix;  0: no, 1: yes */
    OracleMethod,	/* 0: primal, 1: dual */
    OracleRatioTest,	/* 0: standard, 1: Harris */
    OraclePricing,	/* 0: standard, 1: steepest edge */
    Direction,		/* 0: minimize, 1: maximize, set by the Oracle */
    ARGm,		/* -m[0-3] option */
    ARGy,		/* -y or -y- */
    ARGm_set,ARGy_set,ARGp_set,ARGr_set,ARGk_set;
			/* whether these options were used */

  int			/* integer parameters */
    ShowProgress,	/* progress frequency in seconds, 0 means no report */
    RecalculateFacets,	/* after that many iterations do it */
    CheckConsistency,	/* after that many iterations do it */
    OracleOutFreq,	/* oracle output frequency in seconds, >=5 */
    OracleItLimit,	/* iteration limit, >=1000; =0: unlimited */
    OracleTimeLimit,	/* time limit in seconds, >=5; =0: unlimited */
    ProblemColumns,	/* problem columns, set by the Oracle */
    ProblemRows,	/* problem rows, set by the Oracle */
    ProblemObjects,	/* problem objects (dimension), set by the Oracle */
    ARGp,		/* -p T */
    ARGr,		/* -r N */
    ARGk;		/* -k N */

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
    *ProblemName,	/* the problem name, typically the base of the vlp file */
    *ConfigFile,	/* configuration file name */
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
*/
int process_parameters(int argc, const char *argv[]);

/* EOF */

