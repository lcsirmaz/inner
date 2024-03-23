/** glp_oracle.h  -- vertex separation oracle **/

/***********************************************************************
 * This code is part of INNER, a linear multiobjective problem solver.
 *
 * Copyright (C) 2016-2024 Laszlo Csirmaz, https://github.com/lcsirmaz/inner
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
* struct OracleData
*    contains the question to be asked from the oracle and the answer
*    returned. It is allocated by load_vlp(), and never released.
*
* int load_vlp()
*    read the polytope description from a vlp file, store problem
*    dimensions in PARAMS(), create the glpk LP problem instance.
*    Return values:
*      0: file read, memory allocated, glpk LP object initialized
*      1: some error (syntax error, out of bound values, no memory);
*         errors are reported as R_fatal. The vlp file might be left
*         open. The program should abort, no way to recover
*/

typedef struct {
    double *ofacet;		/* ofacet[0..objs-1] the request */
    double *overtex;		/* overtex[0..objs-1] the response */
} OracleData_t;

extern OracleData_t OracleData;
int load_vlp(void);

/**********************************************************************
* Oracle manipulation
*
* int initialize_oracle(void)
*    check the consistency of the vlp problem; compute an initial
*    vertex to OracleData.overtex. Return values:
*      ORACLE_OK     success
*      ORACLE_UNBND  problem is unbounded in some object direction
*      ORACLE_EMPTY  no feasible soultion
*      ORACLE_FAIL   other error condition (limit reached, solver failed)
*
* int ask_oracle(void)
*    the question and the answer is provided in OracleData. Return:
*      ORACLE_OK     success
*      other         some error (limit reached, error, etc)
*/
#define ORACLE_OK	0
#define ORACLE_UNBND	1	/* the projection is unbounded */
#define ORACLE_EMPTY	2	/* the polytope is empty */
#define ORACLE_FAIL	3	/* the oracle failed */

int initialize_oracle(void);
int ask_oracle(void);

/**********************************************************************
* Get oracle statistics
*
* void get_oracle_stat(&callno, &roundno, &time)
*    return the number of oracle calls, total number of iterations,
*    and the used time in 0.01 seconds.
*/
void get_oracle_stat(int *callno,int *roundno,unsigned long *time);

/* EOF */

