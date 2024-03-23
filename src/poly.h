/** poly.h --  combinatorial part using double description method **/

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

/**********************************************************************
* Constant values
*
* BITMAP_64, BITMAP_32
*    bitmaps use either 64 (default) or 32 bit wide words. Use
*    -DBITMAP32 as compiler argument to use 32 bits in bitmaps.
*
* int MAXIMAL_ALLOWED_DIMENSION = 100
*    the maximal dimension we are willing to handle.
*
* int DD_INITIAL_VERTEXNO	= 124
* int DD_INITIAL_FACETNO	= 4092
*    initial number of vertices and facets.
*
* int DD_VERTEX_ADDBLOCK	= 2
* int DD_FACET_ADDBLOCK		= 64
*    space for vertices and facets are allocated in large blocks.
*    These values multiplied by bitmap word size (64 or 32) determine
*    how many new vertices and facets will be added.
*
* size_t DD_HIGHWATER		= 50M
* size_t DD_LOWWATER		= 10M
*    when a block has DD_WIGHWATER more bytes allocated than required,
*    it is reduced back to an excess of DD_LOWWATER
*/

/** maximal dimension we are willing to handle **/
#define MAXIMAL_ALLOWED_DIMENSION	100

/** initial number of vertices/facets. **/
#define DD_INITIAL_VERTEXNO	124
#define DD_INITIAL_FACETNO	4092

/** controlling too large memory blocks, in bytes **/
#define DD_HIGHWATER	((size_t)50e+6)   // 50M
#define DD_LOWWATER	((size_t)10e+6)   // 10M

/** asking space for 128 vertices and 4096 facets **/
#ifdef BITMAP_32		/* 32 bit bitmap blocks */
#define DD_VERTEX_ADDBLOCK	4
#define DD_FACET_ADDBLOCK	128
#else				/* 64 bit bitmaps blocks */
#define DD_VERTEX_ADDBLOCK	2
#define DD_FACET_ADDBLOCK	64
#endif

#ifdef USETHREADS
/************************************************************************
* Threads
*
* int create_threads(void)
*    initialize all threads. Return non-zero in case of an error.
* void stop_threads(void)
*    clear up threads
*/
int create_threads(void);
void stop_threads(void);

#endif /* USETHREADS */

/************************************************************************
* Statistics
*
* Several values are collected during the run of the DD algorithm,
*   which will be printed out when requested. Error conditions, such
*   as inconsistency, or out of memory, are also indicated here.
*
* void get_dd_facetno(void)
*   compute the number of facets of the most recent approximation as
*   well as how many of them is known to be a facet of the final
*   polyhedron. Keeping these values up-to-date at each iteration
*   would be too expensive, so we get them when needed.
*/
typedef struct {
/** statistics **/
int iterations;		    /* number of iterations which is the same as
			       the number of vertices added */
int vertexno;               /* vertices generated so far */
int facetenquires;	    /* number of times a facet was requested */
int probevertex;	    /* number of time vertices were probed */
int vertices_allocated_no;  /* number of times vertex space was extended */
int vertices_allocated;	    /* total number of vertices allocated */
int facets_allocated_no;    /* number of times facet space was extended */
int facets_allocated;	    /* total number of facets allocated */
int facets_compressed_no;   /* times facet compression is called */
int facet_pos;              /* last number of positive facets */
int facet_zero;             /* facets adjacent to this vertex */
int facet_neg;              /* negative facets to be dropped */
int facet_new;              /* new facets added */
int max_facets;		    /* maximal number of intermediate facets */
int max_facetsadded;	    /* maximal number of facets added in an iteration */
double avg_facetsadded;	    /* average number of facets added */
double max_tests;	    /* maximal number of ridge tests by iteration */
double avg_tests;	    /* average number of ridge tests by iterations */
int living_facets_no;	    /* filled by get_dd_facetno() */
int final_facets_no;	    /* filled by get_dd_facetno() */
int memory_allocated_no;    /* number of times memory expanded */
size_t total_memory;	    /* total memory (in bytes) actually allocated */
size_t max_memory;          /* maximum memory allocated so far */
/** warning **/
int instability_warning;    /* number of warnings when recalculating facet eqs */
/** error conditions **/
int numerical_error;	    /* numerical error, data is inconsistent */
int out_of_memory;	    /* out of memory, cannot continue */
int data_is_consistent;     /* in case of memory shortage, indicate data consistency */
} DD_STATS;

extern DD_STATS dd_stats;
void get_dd_facetno(void);

/************************************************************************
* Iterations of the double description algorithm
*
* int init_dd_structure(int maxvertex, int maxfacet)
*    initialize DD algorithm structure allocating space to accomodate
*    the given number of vertices and facets.
*    Return value:
*     0:  initialization is succesful
*     1:  either dimension is too large, or out of memory.
*
* void init_dd(double v[0:dim-1])
*    Add the very first vertex, and create the first inner polytope
*    connecting the vertex with the positive ideal points.
*
* int initial_vertex(double coords[0..dim-1])
* int initial_facet(int final, double coords[0..dim])
*    add initial vertex and facet with the given coords. For a facet,
*    flag final is set if the facet is final. Check that there is enough
*    room for the new item; compute and set incidence. No vertex should
*    be on the negative side of any facet.
*    Return value:
*     0:  OK
*     1:  some error, error message issued.
*
* int get_next_facet(int from, double to[0:DIM])
*    Get the index of a facet of the actual approximation which is
*    not marked as facet of the final polytope. When 'from' is non-
*    negative, return only facets which have index >= from. If
*    from==-1, then return the smallest index, or a random index
*    depending on the parameter RandomFacet. Return -1 if no facet
*    is found.
*    
* void mark_facet_as_final(fno)
*    Mark facet with index fno as facet of the final polytope.
*
* void clear_facet_as_living(fno)
*    This is not a facet in the final polytope
*
* void add_new_vertex(double v[0:dim-1])
*    the main routine, the argument is a new vertex which is to be
*    added to the approximation. Perform all computations. When
*    the routine returns, check error conditions in dd_stats.
*/

/** initialize data structures with the first vertex **/
int init_dd_structure(int maxvertex, int maxfacet);
void init_dd(const double *coords);

/** add initial vertex and facet **/
int add_initial_vertex(const double coords[]);
int add_initial_facet(int final,const double coords[]);

/** index of the next living but not final facet **/
int get_next_facet(int from,double *to);

/** mark the facet as final **/
void mark_facet_as_final(int fno);
/** this facet is not living anymore **/
void clear_facet_as_living(int fno);

/** add a new vertex **/
void add_new_vertex(const double *coords);

/************************************************************************
* Retrieving data, checking, recalculating
*
* int get_vertexnum()
* int get_facetnum()
*    Return the number of non-ideal vertices and facets.
*
* void free_adjacency_lists()
*    After interrupt adjacency lists are not used anymore. Release
*    them to save memory.
*
* int store_vertex(double v[0:dim-1])
*    Called after an interrupt. Check if the argument is among the
*    vertices stored, and add if it was not there. Check memory
*    overflow after return.
*
* int probe_vertex(double v[0:dim-1])
*    Return a score; the vertex with the highest score will be added
*    to the approximation.
*
* void recalculate_facets(void)
*    Go over all facets and recalculate their equations from the list
*    of vertices it is adjacent to. When the routine returns, error
*    conditions should be checked.
*
* int check_bitmap_consistency(void)
*    Check that all facets and vertices have at least DIM adjacent 
*    elements.
*
* int check_consistency(void)
*    Check consistency of the data structure. Return zero if all is
*    fine, otherwise the number of problems detected. 
*/

/** number of vertices and facets **/
int get_vertexnum(void); int get_facetnum(void);

/** release vertex and facet adjacency lists **/
void free_adjacency_lists(void);

/** store a vertex; tell if this is a new one **/
int store_vertex(double *coords);

/** return vertex score **/
int probe_vertex(double *coords);

/** recalculate facet equations **/
void recalculate_facets(void);

/** checking data structure consistency **/
int check_bitmap_consistency(void);
int check_consistency(void);

/************************************************************************
* Report the facets and vertices of the solution
*
* void print_vertex(report_type channel, double coords[])
*    Report the coordinates of a vertex on the given channel. Use
*    fractional format if VertexAsFraction is set, otherwise a floating
*    number. Coordinates are separated by a space, and a terminating
*    newline is added at the end.
*
* void print_vertices(report_type channel)
*    Report all vertices on the given channel. The coordinates of each
*    vertex is printed on a separate line. The first character is a V,
*    the coordinates are printed using print_vertex().
*
* void print_facet(report_type channel, facetno)
*    Report the facet equation, starting with F for final, f for other
*    facets. The equation is scaled first making the coeffs close to
*    integers.
*
* void print_facets(report_type channel)
*    Report the equation of all living facets of on the given channel
*    using print_facet()
*
* void report_memory_usage(report_type ch, int force, char *prompt)
*    report memory usage of the inner approximation algorithm.
*
* void make_checkpoint(void)
*    create the next checkpoint file. Do not report any problem.
*
* void make_dump(void)
*    create a dump of vertices and facets when requested by a signal.
*/

/** reporting vertices **/
void print_vertex(report_type channel, const double coords[]);
void print_vertices(report_type channel); 

/** reporting facets **/
void print_facet(report_type channel,int facetno);
void print_facets(report_type channel);

/** report memory usage **/
void report_memory_usage(report_type channel, int force, const char *prompt);

/* create checkpoint files **/
void make_checkpoint(void);

/* create dump */
void make_dump(void);

/* EOF */


