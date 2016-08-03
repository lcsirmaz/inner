/** report.c  --  report info, error, result **/

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

#include "report.h"
#include "params.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

/***********************************************************************
* Report routines
*
* int check_outfiles(void)
*    check files specified in SaveFile, SaveVertexFile, SaveFacetFile
*    if they are writable. Return 0 if yes, 1 if no.
*
* int checkfile(char *fname)
*    open file 'fname' for writing, and close immediately. Return 0 if
*    successful, 0 if not.
*
* void flush_report(void)
*    if there is a pending output on stdout, flush it.
*
* void report(report_type, format, ...)
*    send the value to the given channel; suppress message depending
*    on PARAMS(). Flush immediately for R_warn and R_info; save on
*    other types. Channels R_savefacet and R_savevertex open files
*    specified in PARAMS().
*
* void close_savefiles(void)
*    close files corresponding to R_savevertex and R_savefacet
*/

static int checkfile(const char *fname)
{FILE *f;
    if(!fname) return 0; /* yes */
    f=fopen(fname,"w"); /* truncate or create */
    if(f){ fclose(f); return 0; }
    report(R_fatal,"Error: cannot open output file %s for writing.\n",fname);
    return 1;
}

int check_outfiles(void)
{
    return checkfile(PARAMS(SaveFile)) ||
           checkfile(PARAMS(SaveVertexFile)) ||
           checkfile(PARAMS(SaveFacetFile));
}

/***********************************************************************
*
* int pending_output
*    non-zero if there is a pending output on stdout
*
* FILE *savefile, *savevertexfile, *savefacetfile, *chkfile
*    either NULL or the opened stream handle for that channel.
*
* int open_vertexfile(void)
*    if savevertexfile is not open, open it. Take care if it is
*    the same as savefacetfile, and that has been opened.
*
* int open_facetfile(void)
*    if savefacetfile is not open, open it. Take care if it is
*    the save as savevertexfile, and that has been opened.
*/
static int pending_output=0;
static FILE *savefile=NULL, *savevertexfile=NULL, 
            *savefacetfile=NULL, *chkfile=NULL;

static int open_vertexfile(void)
{   if(savevertexfile) return 1;
    if(PARAMS(SaveFacetFile)
       && strcmp(PARAMS(SaveFacetFile),PARAMS(SaveVertexFile))==0
       && savefacetfile){
         savevertexfile=savefacetfile; return 1;
    }
    savevertexfile=fopen(PARAMS(SaveVertexFile),"w");
    return savevertexfile!=NULL;
}

static int open_facetfile(void)
{   if(savefacetfile) return 1;
    if(PARAMS(SaveVertexFile)
       && strcmp(PARAMS(SaveVertexFile),PARAMS(SaveFacetFile))==0
       && savevertexfile){
        savefacetfile=savevertexfile; return 1;
    }
    savefacetfile=fopen(PARAMS(SaveFacetFile),"w");
    return savefacetfile!=NULL;
}

/** report(report_type, fmt, ...) **/
void report(report_type level, const char *fmt,...)
{va_list arg; int y,flush;
    y=0; flush=0;
    switch(level){
      case R_info: if(PARAMS(MessageLevel)>2){ y=1; flush=1; } break;
      case R_fatal:
      case R_txt:  y=1; break;
      case R_err:  if(PARAMS(MessageLevel)>0) y=1; break;
      case R_warn: if(PARAMS(MessageLevel)>1) { y=1; flush=1; } break;
      default:     y=-1; break;
    }
    if(y>=0){
        if(y){
            va_start(arg,fmt); vprintf(fmt,arg); va_end(arg);
            pending_output=1;
            if(flush){ pending_output=0; fflush(stdout); }
        }
        return;
    }
    if(level==R_savefacet){
        if(PARAMS(SaveFile) && PARAMS(SaveFacets)){
            if(!savefile) savefile=fopen(PARAMS(SaveFile),"w");
            if(savefile){
              va_start(arg,fmt); vfprintf(savefile,fmt,arg); va_end(arg);
            }
        }
        if(PARAMS(SaveFacetFile) && open_facetfile()){
            va_start(arg,fmt); vfprintf(savefacetfile,fmt,arg); va_end(arg);
        }
    } else if(level==R_savevertex){
        if(PARAMS(SaveFile) && PARAMS(SaveVertices)){
            if(!savefile) savefile=fopen(PARAMS(SaveFile),"w");
            if(savefile){
              va_start(arg,fmt); vfprintf(savefile,fmt,arg); va_end(arg);
            }
        }
        if(PARAMS(SaveVertexFile) && open_vertexfile()){
            va_start(arg,fmt); vfprintf(savevertexfile,fmt,arg); va_end(arg);
        }
    } else if(level==R_chk){
        if(chkfile){
          va_start(arg,fmt); vfprintf(chkfile,fmt,arg); va_end(arg);
        }
    }
}

void close_savefiles(void)
{   if(savefile) { fclose(savefile); savefile=NULL; }
    if(savevertexfile){
        if(savevertexfile==savefacetfile) savefacetfile=NULL;
        fclose(savevertexfile); savevertexfile=NULL;
    }
    if(savefacetfile){ fclose(savefacetfile); savefacetfile=NULL; }
}

void flush_report(void)
{   if(pending_output){ pending_output=0; fflush(stdout); } }

void open_checkpoint(int version)
{static char *buf=NULL;
    if(version<0 || version>999 ) version=0;
    if(chkfile){ fclose(chkfile); chkfile=NULL;}
    if(buf==NULL){
        buf=malloc(10+strlen(PARAMS(CheckPointStub)));
        if(!buf) return; // out of memory
    }
    sprintf(buf,"%s%03d.chk",PARAMS(CheckPointStub),version);
    chkfile=fopen(buf,"w"); // truncate or create
}

void close_checkpoint(void)
{   if(chkfile){ fclose(chkfile); chkfile=NULL; } }

/* EOF */

