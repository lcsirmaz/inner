/** main.c -- the matin program **/

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

#include "main.h"
#include "params.h"
#include "inner.h"
#include "report.h"

/*******************************************************************
* signal handling 
*
* volatile int dobreak
*    inticates how many timer SIGNAL was sent
* void set_ctrl_c()
*    whenever SIGNAL is received, increase the value of dobreak.
*    This is to be checked in the main loops.
*/

#include <signal.h>
#include <stdio.h>

volatile int dobreak=0;
static void siginthandler(int signo)
{   if(signo==INNER_SIGNAL) dobreak++; }

static int set_signal(void)
{struct sigaction act;
    sigemptyset(&act.sa_mask); // which signals to be blocked
    act.sa_handler = siginthandler;
    act.sa_flags = 0; // do not reset to SIGN_DFL
    if(sigaction(INNER_SIGNAL, &act, NULL)){
        report(R_fatal,"Cannot set " mkstringof(INNER_SIGNAL) " signal, aborting\n");
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
    if(set_signal()) return 1; // handle intterupts
    switch(inner()){ // execute the algorithm
      case 0:  return 0; /* job done */
      case 1:  return 1; /* data error before algorithm started */
      case 2:  return 2; /* error in the algorithm */
      default: break;    /* interrupted, OK, error, aborted */
    }
    return 3;
}

/* EOF */

