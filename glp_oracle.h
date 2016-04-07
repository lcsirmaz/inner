/** glp_oracle.h  -- vertex separation oracle **/

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
*  The *vertex separation oracle* has a (hidden) polytope. A question is
*  a vector 'h' (a hyperplane cutting into the polytope). The answer is a
*  vertex 'v' of the polytope which minimizes the inner product 'v*h'
*  This routine realizes this oracle using a patched version of the Gnu
*  Linear Program Kit (glpk). The hidden polytope is defined in a vlp
*  file.
*/

/**********************************************************************
* 
* VertexOracleData: struct containing problem dimensions, the question
*    to be asked from the oracle, and the answer it returns. The structure
*    is allocated by read_vlp(), and never released.
*/

typedef struct {
    double *ofacet;		/* ofacet[0..objs-1] the request */
    double *overtex;		/* overtex[0..objs-1] the response */
} VertexOracle_t;

extern VertexOracle_t VertexOracleData;

/**********************************************************************
* read_vlp: read the polytope description from a vlp file. Store 
*   dimensions in VertexOracleData. Create the LP problem instance.
* int read_vlp(void)
*  Return value:
*    0: file read, memory allocated, glpk LP object initialized
*    1: some error (syntax error, out of bound values, no memory);
*      errors are reported as R_fatal. The vlp file might be left
*      open. The program should abort, no way to recover.
*
* set_oracle_parameters: set the parameters of the LP solver from
*   the configuration (primal/dual, time bound, verbosity, etc)
* void set_oracle_parameters(void)
*
*/
int read_vlp(void);

void set_oracle_parameters(void);

/**********************************************************************
* Ask the oracle
*
* int ask_oracle(void)
*    the question and the answer is provided in VertexOracleData.
*
* Return values:
*   ORACLE_OK     the minimal vertex is stored in VertexOracleData,
*                 coordinates are rounded to the nearest rational value
*                 when "RoundVertices" is set.
*   ORACLE_UNBND  the polytope is not bounded from below
*   ORACLE_EMPTY  the polytope is empty (no feasible solution)
*   ORACLE_LIMIT  either time or iteration limit is reached
*   ORACLE_FAIL   the LP solver failed to solve the problem
*/
#define ORACLE_OK	0
#define ORACLE_UNBND	1	/* the projection is unbounded */
#define ORACLE_EMPTY	2	/* the polytope is empty */
#define ORACLE_LIMIT	3	/* iteration or time limit reached */
#define ORACLE_FAIL	4	/* the oracle failed */

int ask_oracle(void);

/**********************************************************************
* Get LP solver iterations
*
* int get_oracle_rounds()
*    how many iterations (rounds) the LP solver did to solve this
*    particular problem. 
*/
int get_oracle_rounds(void);

/* EOF */

