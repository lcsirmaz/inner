/** version.h -- version and copyright information **/

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

/* Version and copyright */
#define VERSION_MAJOR	2
#define VERSION_MINOR	15
#ifdef USETHREADS
#define VERSION_STRING	"threaded version " mkstringof(VERSION_MAJOR.VERSION_MINOR) "T"
#else
#define VERSION_STRING	"version " mkstringof(VERSION_MAJOR.VERSION_MINOR)
#endif

#define COPYRIGHT	\
"Copyright (C) 2016-2024 Laszlo Csirmaz, https://github.com/lcsirmaz/inner"

/* EOF */

