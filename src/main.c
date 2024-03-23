/** main.c -- the main program **/

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

#include "main.h"
#include "params.h"
#include "inner.h"
#include "report.h"

/*******************************************************************
* signal handling 
*
* volatile int dobreak, dodump
*    inticates how many timer a SIGNAL was sent
* void set_signals()
*    whenever INNER_SIGNAL or DUMP_SIGNAL is received, increase
*    the value of the above values. They are checked in the main
*    loop of the program.
*/

#include <signal.h>
#include <stdio.h>

volatile int dobreak=0, dodump=0;
static void siginthandler(int signo)
{   if(signo==INNER_SIGNAL) dobreak++;
    if(signo==DUMP_SIGNAL)  dodump++;
}

inline static int set_signals(void)
{struct sigaction act;
    sigemptyset(&act.sa_mask); // which signals to be blocked
    act.sa_handler = siginthandler;
    act.sa_flags = 0; // do not reset to SIGN_DFL
    if(sigaction(INNER_SIGNAL, &act, NULL)){
        report(R_fatal,"Cannot set signal " mkstringof(INNER_SIGNAL) ", aborting\n");
        return 1;
    }
    if(sigaction(DUMP_SIGNAL, &act, NULL)){
        report(R_fatal,"Cannot set signal " mkstringof(DUMP_SIGNAL) ", aborting\n");
        return 1;
    }
    return 0;
}

/*******************************************************************
* the main program
*     read parameters, input file
*     execute the inner approximation algorithm
*/

int main(int argc, const char *argv[])
{int r;
    r=process_parameters(argc,argv); // read parameters
    if(r==1) return 0; /* job done */
    if(r<0) return 1;  /* data error */
    if(set_signals()) return 1; // handle signals
    switch(inner()){ // execute the algorithm
      case 0:  return 0; /* job done */
      case 1:  return 1; /* data error before algorithm started */
      case 2:  return 2; /* no solution */
      case 3:  return 3; /* problem unbounded */
      case 4:  return 4; /* error in the algorithm */
      default: break;    /* interrupted, OK, error, aborted */
    }
    return 5;
}

/* EOF */

