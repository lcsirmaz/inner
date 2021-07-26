/** poly.h --  combinatorial part using double description method **/

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

/**********************************************************************
* Constant values
*
* BITMAP_64, BITMAP_32
*    bitmaps use either 64 (default) or 32 bit wide words. Use
*    -DBITMAP32 as compiler argument to use 32 bits in bitmaps.
*
* int MAXIMAL_ALLOWED_DIMENSION
*    the maximal dimension we are willing to handle.
*
* int DD_INITIAL_VERTEXNO
* int DD_INITIAL_FACETNO
*    initial number of vertices and facets.
*
* int DD_VERTEX_ADDBLOCK
* int DD_FACET_ADDBLOCK
*    space for vertices and facets are allocated in large blocks.
*    These values multiplied by bitmap word size (64 or 32) determine
*    how many new vertices and facets will be added.
*/

/** maximal dimension we are willing to handle **/
#define MAXIMAL_ALLOWED_DIMENSION	100

/** initial number of vertices/facets. **/

#define DD_INITIAL_VERTEXNO	124
#define DD_INITIAL_FACETNO	4092

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
*
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
*   computes the number of facets of the most recent approximation as
*   well as how many of them is known to be a facet of the final
*   polyhedron. Keeping these values up-to-date at each iteration
*   would be too expensive.
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
int facet_new;              /* facets added */
int max_facets;		    /* maximal number of intermediate facets */
int max_facetsadded;	    /* maximal number of facets added in an iteration */
double avg_facetsadded;	    /* average number of facets added */
double max_tests;	    /* maximal number of ridge tests by iteration */
double avg_tests;	    /* average number of ridge tests by iterations */
int living_facets_no;	    /* filled by get_dd_facetno() */
int final_facets_no;	    /* filled by get_dd_facetno() */
int memory_allocated_no;    /* number of times memory expanded */
size_t total_memory;	    /* total memory (in bytes) allocated */
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
* int init_dd_structure(int dim, int maxvertex, int maxfacet)
*    initialize DD algorithm structure allocating space to accomodate
*    the given number of vertices and facets.
*    Return value:
*     0:  initialization is succesful
*     1:  either dimension is too large, or out of memory.
*
* int init_dd(int dim, double v[0:dim-1])
*    initialize the DD algorithm, supplying the very first vertex 
*    allocation the initial number of vertices and facets.
*    The algorithms is tailored to MOLP minimization; the positive
*    endpoints of coordinate axes (the *ideal* vertices) are assumed
*    to be feasible solutions. Arguments:
*      dim:         the dimension of the space.
*      v[0:dim-1]:  coordinates of the very first extremal solution.
*    Return value:
*     0:  initialization is successful
*     1:  either the dimension is too large, or out of memory.
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
* int get_next_facet(int from)
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
* void get_facet_into(int fno, double eq[0:dim])
*    Retrieve the equation of the facet with index fno. In the
*    internal representation all coefficients in eq[0:dim-1] are
*    non-negative and sum to 1.0. When returning the equation it
*    is scaled and after scaling values closer to and integer than
*    ScaleEps are rounded up to that integer.
*
* void add_new_vertex(double v[0:dim-1])
*    Specify the vertex which enlarges the inner approximation. When
*    the routine returns, check error conditions in dd_stats.
*/

/** initialize data structures with the first vertex **/
int init_dd_structure(int dim, int maxvertex, int maxfacet);
int init_dd(int dim, double *coords);

/** add initial vertex and facet **/
int initial_vertex(double coords[]);
int initial_facet(int final, double coords[]);

/** release vertex and facet adjacency lists **/
void free_adjacency_lists(void);

/** index of the next living but not final facet **/
int get_next_facet(int from);

/** mark the facet as final **/
void mark_facet_as_final(int facetno);
/** this facet is not living anymore **/
void clear_facet_as_living(int facetno);
/** retrieve facet equation **/
void get_facet_into(int facetno, double *faceteq);

/** add a new vertex **/
void add_new_vertex(double *coords);

/************************************************************************
* Retrieving data, checking, recalculating
*
* int vertex_num()
* int facet_num()
*    Return the number of non-ideal vertices and facets.
*
* int store_vertex(double v[0:dim-1])
*    If the algorithm is stopped (it runs for too long time), 
*    extremal solutions to the MOLP can still be generated by
*    asking the Oracle all facets of the last approximation.
*    the returned vertex is then stored as a final vertex.
*    Return 1 if this is a new vertex not seen before, so it
*    can be reported. Error conditions should be checked
*    after the routine returns.
*
* int probe_vertex(double v[0:dim-1])
*    Return the number of facets which would be thrown away if this 
*    vertex were added to the approximation.
*
* void recalculate_facets(void)
*    Go over all facets and recalculate their equations from the list
*    of vertices it is adjacent to. When the routine returns, error
*    conditions should be checked.
*
* int check_consistency(void)
*    Check consistency of the data structure. Return zero if all is
*    fine, otherwise the number of problems detected. 
*/

/** number of vertices and facets **/
int vertex_num(void); int facet_num(void);

/** store a vertex; tell if this is a new one **/
int store_vertex(double *coords);

/** return the number of negative facets **/
int probe_vertex(double *coords);

/** recalculate facet equations **/
void recalculate_facets(void);

/** checking data structure consistency **/
int check_consistency(void);

/************************************************************************
* Report the facets and vertices of the solution
*
* void print_vertex(report_type channel, double dir, double coords[])
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
* void report_memory_usage(void)
*    report memory usage of the inner approximation algorithm.
*
* void make_checkpoint(void)
*    create the next checkpoint file. Do not report any problem.
*
* void make_dump(void)
*    create a dump of vertices and facets when requested by a signal.
*/

/** reporting vertices **/
void print_vertex(report_type channel, double dir, const double coords[]);
void print_vertices(report_type channel); 

/** reporting facets **/
void print_facet(report_type channel,int facetno);
void print_facets(report_type channel);

/** report memory usage **/
void report_memory_usage(void);

/* create checkpoint files **/
void make_checkpoint(void);

/* create dump */
void make_dump(void);

/* EOF */


