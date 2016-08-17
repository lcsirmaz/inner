/** version.h -- version and copyright information **/

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

/* Version and copyright */
#define VERSION_MAJOR	2
#define VERSION_MINOR	6
#ifdef USETHREADS
#define VERSION_STRING	"T" mkstringof(VERSION_MAJOR.VERSION_MINOR)
#else
#define VERSION_STRING	mkstringof(VERSION_MAJOR.VERSION_MINOR)
#endif

#define COPYRIGHT	\
"Copyright (C) 2016 Laszlo Csirmaz, Central European University, Budapest"

/* EOF */

