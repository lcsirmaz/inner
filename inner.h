/** inner.h  -- bounds and default values **/

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
* the inner approximation algorithm
*
* int inner(void)
*    execute the inner approximation algorithm for the prepared
*    data.
*    Return value:
*      0   program terminated normally
*      1   data error before the algorithm started
*      2   error during algorithm execution
*      3   interrupted, postprocessing terminated normally
*      4   interrupted, error in postprocessing
*      5   interrupted, postprocessing aborted
*/

int inner(void);

/* EOF */

