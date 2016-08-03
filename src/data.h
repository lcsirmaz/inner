/* data.h -- parse and read data files */

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
* int init_reading(char *fname)
*   start reading data from the given file. Return values:
*     1: some error, error message issued (out of memory, cannot open file)
*     0: OK
*
* int nextline( int *type)
*   read the next line into the internal buffer. Return values:
*     1: next line read, *type contains the type of the line:
*         1: V, 2: F, 3: f, 4: N
*     0: no more lines, file handle is closed
*
* int parseline(int num, double to[0..num-1])
*   parse the content of the last line. It should contain exactly num
*   double values which are stored to the given place. Return values:
*      1: some error (no message is given)
*      0: all values are read and stored
*/

int init_reading(const char *fname);
int nextline(int *type);
int parseline(int num, double *to);

/* EOF */


