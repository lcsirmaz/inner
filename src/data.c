/** data.c parse and read data files **/

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "report.h"
#include "data.h"

/**********************************************************************
* Read vertices and facets from a file
*
* int LINELEN = 1000
*   maximum length of a line
* char line[LINELEN]
*   next line from the file
* FILE *handle
*   filehandle, open while the file is not over
* char *filename
*   the original file name to show in error messages
*/
#define LINELEN		1000
static char *line=NULL;
static FILE *handle=NULL; /* filehandle */
static char *filename=NULL;

/* int nextline(int *type)
*   read the next line starting with 'V', 'F', or 'f' into line[]. If
*   no more lines, close the filehandle.
*   Return values:
*     1: next line is in line[]; *type is set to the type of the line
*     0: EOF, filehandle is closed.
*/
int nextline(int *type)
{int i,sp,ch;
    if(handle==NULL) return 0;
    i=0;sp=0; memset(line,0,LINELEN+1);
    while((ch=getc(handle))>=0){
        if(ch=='\n'){ if(i<=0) {sp=0;i=0; continue; } return 1; }
        if(i<0) continue; // skip this line
        if(i==0 && sp==0){
            if(ch=='V'){ *type=1; sp=1; }
            else if(ch=='F'){ *type=2; sp=1; }
            else if(ch=='f'){ *type=3; sp=1; }
            else if(ch=='N'){ *type=4; sp=1; }
            else {i=-1; }
            continue;
        }
        if(ch==' '||ch=='\t'){ sp=1; continue; }
        if(ch<=0x20 || ch>126) continue; // ignore these characters
        if(sp && i>0 && i<LINELEN){line[i]=' '; i++; }
        sp=0; if(i<LINELEN){line[i]=ch; i++; }
    } /* EOF */
    if(i>0) return 1;
    fclose(handle); handle=NULL;
    return 0;
}
/* int parse_num(char *str, double *v)
*   check if str[] starts with a (floating) number. If yes, store the
*   value in *v, and return the length of the matching part.
*   Return value:
*     0: no more value, end of string hit
*    -1: syntax error
*    >0: the length of the matching part
*/
static int parse_num(char *str, double *v)
{int i; char *s;
    if(!*str) return 0; // no more input
    for(i=0,s=str;('0'<=*s && *s<='9')|| *s=='+'||*s=='-'
                ||*s=='.'||*s=='e'; i++,s++);
    if(i==0) return -1; // syntax error
    if(sscanf(str,"%lf",v)!=1) return -1; //syntax error
    return i;
}
/* int parse_real(char *str, double *v)
*   check whether str[] starts with a number, or number/number. If yes,
*   store the value in *v and return the length of the matching part.
*   Return value is the same as for the previous routine.
*/
static int parse_real(char *str, double *v)
{int i,ret; double d;
    i=0; if(*str==' '){ str++; i++; }
    ret=parse_num(str,v);
    if(ret<=0) return ret; // 0 or -1
    i+=ret; str+=ret;
    if(*str=='/'){
        str++; i++; d=*v;
        ret=parse_num(str,v);
        if(ret<=0) return -1; // syntax error
        i+=ret; *v = d/(*v);
    }
    return i;
}
/* int parseline(int num, double to[0..num-1])
*   stores the values at the next line to the array provided
*   Return values:
*    0: OK
*    1: some error
*/
int parseline(int num, double *to/*0..num-1*/)
{char *str; int i,ret;
    str=line;
    for(i=0;i<num;i++,to++){
        ret=parse_real(str,to);
        if(ret<=0){// error
            report(R_fatal,"Error parsing data line in file %s\n%s\n",filename,line);
            return 1;
        }
        str+=ret;
    }
    if(*str){ // extra characters at the end
        report(R_fatal,"Extra characters at the end of data line in file %s\n%s\n",filename,line);
        return 1;
    }
    return 0;
}
/* int init_reading(char *fname)
*   start reading data from the given file. Return values:
*     1: some error, error message issued (out of memory, cannot open file)
*     0: OK
*/
int init_reading(const char *fname)
{
    if(!line) line=malloc(LINELEN+4);
    if(!line){ report(R_fatal,"Out of memory (init_reading)\n"); return 1; }
    handle=fopen(fname,"r");
    if(!handle){ report(R_fatal,"Cannot open file %s for reading\n",fname); return 1; }
    filename=strdup(fname); // for error messages
    return 0;
}

/* EOF */


