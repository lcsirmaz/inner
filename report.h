/** report.h  -- report info, error, results **/

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
* Reporting
*
* report_type
*    different type of reports: info, warn, err, fatal
*
* void report(report_type, format, ...)
*    depending on report type and parameter setting, send or suppress
*    message to stdout or to a result file.
*
* void flush_report(void)
*    flush pending messages at stdout
*
* void close_savefiles(void)
*    close files corresponding to R_savevertex and R_savefacet
*
* int check_outfiles(void)
*    check that all result files are writable. Do it before starting
*    any serious computation as we don't want "cannot write file"
*    after a day's work.
*/

/* report type */
typedef enum {
R_fatal,	/* fatal error */
R_err,		/* error */
R_warn,		/* warning  + flush */
R_txt,		/* report message */
R_info,		/* info + flush */
R_savefacet,	/* save facets, go to result file */
R_savevertex	/* save vertices, go to result file */
} report_type;

/** report the message to the given channel */
void report(report_type level, const char *fmt,...)
#ifdef __GNUC__
__attribute__((format(printf,2,3)))
#endif
;

/** close result files **/
void close_savefiles(void);

/** flush stdout when there is a pending message **/
void flush_report(void);

/** check files if they are writable **/
int check_outfiles(void);

/* EOF */

