/** poly.c  --  combinatorial part using double description method **/

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

/*
* This part of INNER is the combinatorial vertex enumeration algorithm.
* We maintain an inner approximation of the final polytope. In each stage
* a new facet of the approximation is chosen, and the oracle is asked
* whether it is a facet of the final polytope. If yes, it is marked as
* "final". If not, the oracle returns a vertex of the final polytope 
* on the facet's negative side, and the approximation is enlarged.
*
* The algorithm starts with a d dimensional simplex which has d ideal
* vertices at the positive end of the axes, and a "real" vertex. As
* this simplex is spanned by vertices of the final polytope, each facet
* (final and intermediate) has non-negative normal.
*
* Both the vertices and the facets of the inner approximation is stored
* together with adjacency lists which are handled as bitmaps. Using
* bitmaps allows fast operation on vertices and facets.
*/         

#include <stdio.h>
#include <stdint.h>	/* uint32_t, uint64_t */
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "report.h"
#include "poly.h"
#include "params.h"

/************************************************************************
* Memory management routines
*
* Memory is arranged in slots; a slot consists of several blocks of
* given blocksize. When reallocating, either the number of blocks, the
* blocksize, or both can change. When blocksize changes, the old content
* is adjusted to match the new blocksize.
*
* enum memtype_t
*   the managed memory slots
* struct MEMSLOT memory_slots[]
*   static array containing for each slot the actual blocksize, number
*   of blocks, and a pointer to the actual location. The location can
*   change as a result of "defragmentation".
*
* void report_memory_usage(void)
*   report for each slot the blocksize, number of blocks, total memory
*   used, and the actual pointer.
*
* void allocmem(memslot_t slot, int blockno, int blocksize)
*   allocates the requested memory the the slot.
*
* void requestmem(memslot_t slot, int blockno, int blocksize)
*   records the requested new block count and size.
*   tries to allocate more memory if necessary. Return 1 if cannot
*   get more memory, otherwise return 0.
*
* int reallocmem(void)
*   allocate requested memory; if successful, adjust blocks to the
*   given count and size. Return 1 if out of memory, otherwise 0.
*
* void yalloc(type,slot,n,bsize)
*   request n blocks at the given slot, each block is an array of
*   bsize elements of the given type.
*
* void yrequest(type,slot,n,bsize)
*   request more memory at the given slot.
*/

typedef enum { /* managed memory slots */
M_vertexcoord	= 0,	/* vertex coordinates; StoreVertices */
M_vertexadj,		/* adjacency list of vertices; VertexAdjList*/
M_vertexwork,		/* vertex bitmap; VertexWork */
M_facetcoord,		/* facet coordinates; StoreFacets */
M_facetadj,		/* adjacency list of facets; FacetAdjList */
M_facetliving,		/* facet bitmap; FacetLiving */
M_facetfinal,		/* facet bitmap; FacetFinal */
M_facetpos,		/* facet bitmap, FacetPositive */
M_facetneg,		/* facet bitmap, FacetNegative */
M_facetwork,		/* facet bitmap; FacetWork */
M_facetdist,		/* facet distances; StoreFacetDists */
M_vertexarray,		/* calculating facet equations; VertexArray */
M_MSLOTSTOTAL		/* total number of managed memory */
} memslot_t;

typedef struct { /* memory slot structure */
  size_t blocksize;	/* block size (in bytes) */
  size_t blockno;	/* number of blocks */
  size_t newblocksize;	/* new block size */
  size_t newblockno;	/* new block count */
  size_t rsize;		/* real size in bytes */
  void   *ptr;		/* the actual value */
} MEMSLOT;

static MEMSLOT memory_slots[M_MSLOTSTOTAL]; /* memory slots */
#define OUT_OF_MEMORY	dd_stats.out_of_memory

/** report memory usage **/
void report_memory_usage(void)
{int i; MEMSLOT *ms;
    report(R_txt,"vertex and facet storage%s:\n",
        OUT_OF_MEMORY ? " (out of memory)":"");
    for(i=0,ms=&memory_slots[0];i<M_MSLOTSTOTAL;i++,ms++){
        report(R_txt,"M   slot %2d: blocksize=%6u, no=%6u, total=%9u, ptr=%p\n",
           i,ms->blocksize,ms->blockno,ms->rsize,ms->ptr);
    }
}

/* allocate memory to a given slot
*  allocmem(slot,nno,nsize)
*   slot:  the requested slot
*   nno:   number of blocks
*   nsize: block size
*/
static void allocmem(memslot_t slot, size_t nno, size_t nsize)
{size_t total; MEMSLOT *ms;
    if(OUT_OF_MEMORY) return;
    ms=&memory_slots[slot];
    ms->newblocksize=0;
    ms->newblockno=0;
    total=nno*nsize;
    ms->blocksize=nsize;
    ms->blockno=nno;
    ms->rsize=total;
    dd_stats.total_memory += total;
    ms->ptr=malloc(total);
    if(!ms->ptr){
        report(R_fatal,
           "Allocmem: out of memory for slot=%d, blocksize=%d, n=%u\n",
           slot,nsize,nno);
        OUT_OF_MEMORY=1;
        return;
    }
    memset(ms->ptr,0,total);
    return;
}

/* request more memory for a given slot
*  int requestmem(slot,nno,nsize)
*   slot:  the requested slot
*   nno:   number of blocks
*   nsize: block size
*/
static inline void requestmem(memslot_t slot, size_t nno, size_t nsize)
{MEMSLOT *ms;
    ms=&memory_slots[slot];
    if(ms->blocksize==nsize && nno<=ms->blockno) return;
    ms->newblocksize=nsize;
    ms->newblockno=nno;
}    

/* ask more memory and adjust the block structure
* int reallocmem(void)
*   return 0 if OK, 1 if out of memory
*/
static int reallocmem(void)
{MEMSLOT *ms; int j,success; size_t total; void *ptr;
    for(j=0,success=1,ms=&memory_slots[0];
        success && j<M_MSLOTSTOTAL; j++,ms++) if(ms->newblocksize){
        total=ms->newblocksize*ms->newblockno;
        if(ms->rsize<total){
            dd_stats.memory_allocated_no++;
            ptr=realloc(ms->ptr,total);
            if(ptr){
                ms->ptr=ptr; 
                dd_stats.total_memory += total-ms->rsize;
                ms->rsize=total;
            }
            else { success=0; }
        }
    }
    // one can try to defragment memory by writing out all
    // content, then freeing and reallocating memory
    if(!success){ OUT_OF_MEMORY=1; return 1; }
    // adjust block structure
    for(j=0,ms=&memory_slots[0];
        j<M_MSLOTSTOTAL; j++,ms++) if(ms->newblocksize){
        if(ms->blocksize < ms->newblocksize){ // block size grow
            size_t i,n; char *bo,*bn; // pointer to old a new
            n=ms->newblockno; if(ms->blockno<n) n=ms->blockno;
            bo=((char*)ms->ptr)+(n*ms->blocksize);
            bn=((char*)ms->ptr)+(n*ms->newblocksize);
            for(i=n;i>0;i--){
                bo -= ms->blocksize; bn -= ms->newblocksize;
                if(i>1) memmove(bn,bo,ms->blocksize);
                // clear the rest
                memset(bn+ms->blocksize,0,ms->newblocksize-ms->blocksize);
            }
        } else if(ms->blocksize > ms->newblocksize) { // block size shrank
            size_t i,n; char *bo,*bn; // pointers to old and new
            n=ms->newblockno; if(ms->blockno<n)n=ms->blockno;
            bo=(char*)ms->ptr; bn=(char*)ms->ptr;
            for(i=0;i<n;i++){
                if(i) memmove(bn,bo,ms->newblocksize);
                bo += ms->blocksize;
                bn += ms->newblocksize;
            }
        }
        ms->blocksize=ms->newblocksize;
        if(ms->blockno < ms->newblockno){ // clear the rest
            memset(((char*)ms->ptr)+ms->blockno*ms->blocksize,0,
                 (ms->newblockno-ms->blockno)*ms->blocksize);
        }
        ms->blockno=ms->newblockno;
        ms->newblocksize=0;
        ms->newblockno=0;
    }
    return 0;
}

#define yalloc(type,slot,n,bsize)	\
    allocmem(slot,n,(bsize)*sizeof(type))
#define yrequest(type,slot,n,bsize)	\
    requestmem(slot,n,(bsize)*sizeof(type))

/***********************************************************************
* The combinatorial Double Description method
*
* The algorithm works in DIM dimensions, and maintains the vertices and
* facets of the inner approximation (convex, nonempty) polytope. 
*
* vertices of the polytope:
*   The first DIM vertices are *ideal* vertices in the positive endpoints
*   of the coordinate directions. All other vertices are *standard* (or
*   real) vertices. Each vertex has
*   o  DIM coordinates, and
*   o  the adjacency list of facets stored as a bitmap.
*
* facets of the polytope:
*   The first facet is the *ideal* facet adjacent to all ideal vertices.
*   Each facet has
*   o  DIM+1 coordinates. The first DIM coordinates has value >=0 and
*      these values add up to 1.0 -- this normalizes the facet equation.
*      The only exception is the ideal facet, where the first DIM 
*      coordinates are 0.0, and the last coordinate is 1.0.
*   o  the adjacency list vertices on that facet stored as a bitmap
*   o  for each facet its distance from the newly added vertex v; this
*      is computed as
*         distance = f[1]*v[i]+f[2]*v[2]+...+f[DIM]*v{DIM]+f[DIM+1];
*      the polytope is on the >=0 size of each facet.
*   o  two attributes (stored as bitmaps) telling whether the facet is
*      a living one (has not been deleted earlier), and whether this
*      facet is known to be the facet of the final polytope.
*
* All indices start at 0, including vertices, facets, coordinates.
*/

static int DIM; // problem dimension = number of objectives

/**********************************************************************
* B I T M A P S
*
* A *bitmap* is an array of 64 bit (or 32 bit) unsigned integers
* accommodating the requested number of bits. Bitmaps are used as 
* adjacency lists, and also to store vertex and facet attributes.
*
* BITMAP_t, BITMAP0, BITMAP1, packshift, packmask
*   Bitmaps are declared as type BITMAP_t *B. BITMAP0 and BITMAP1
*   are zero and 1 in BITMAP_t type. To get/set "index" bit in
*   bitmap B use packshift and packmask as follows:
*     (B[index>>packshift]>>(index&packmask)) & 1
*     B[index>>packshift] |= BITMAP1 << (index&packmask)
*/

#ifdef BITMAP_32		/* 32 bit bitmap */
typedef uint32_t BITMAP_t;
#define packsizelog	2	/* sizeof(BITMAP_t)== 1<<packsizelog */

#else				/* 64 bit bitmap */
typedef uint64_t BITMAP_t;
typedef union {BITMAP_t u; struct { uint32_t i1,i2;}i; } SPLIT_BITMAP_t;
#define packsizelog	3	/* sizeof(BITMAP_t)== 1<<packsizelog */
#endif

#define packsize	(1<<packsizelog)
#define packshift	(packsizelog+3) 	/* in bits */
#define packmask	((1<<packshift)-1)
#define BITMAP1		((BITMAP_t)1u)    	/* BITMAP_t unsigned 1 */
#define BITMAP0		((BITMAP_t)0u)		/* BITMAP_t unsigned zero */

/********************************************************************
* V E R T I C E S
*
* int VertexSize
*    number of double values to store the coordinates of a vertex
* int NextVertex
*    the index of the first free slot to store a new vertex. The actual
*    number of vertices is NextVertex, and they are indexed from 0 to
*    NextVertex-1
* int MaxVertices
*    the maximal available slot for vertices; NextVertex <= MaxVertices.
*    Vertex bitmaps can store up to MaxVertices bits.
* int VertexBitmapBlockSize
*    number of BITMAP_t words in a vertex bitmap
*
* double *StoreVertices
*    where vertices are stored; each vertex is a block of VertexSize
*    doubles, and there is a room for MaxVertices vertex
* BITMAP_t *FacetAdjList, *VertexWork
*    bitmaps with VertexBitmapBlockSize blocksize. FacetAdjList stores 
*    the adjacency list of facets in MaxFacets blocks. VertexWork is a
*    vertex attribute with a single block storing MaxVertices bits.
* double *VertexArray
*    store (the homogeneous coordinates) of vertices adjacent to a 
*    facet; it is used to compute the facet equation
*
* double *VERTEX_coord(vno)
*    The block of VertexSize coordinates of the vertex vno
* BITMAP_t *VERTEX_adj(vno)
*    the block of the facet adjacency list of vertex vno
*
*********************************************************************
* F A C E T S
*
* int FacetSize
*    number of double values to store a facet equation
* int NextFacet
*    the index of the fist empty slot for a new facet
* int MaxFacets
*    number of bits in a facet bitmap
* int AheadFacets
*    maximal number of slots storing a facet
* int FacetBitmapBlockSize
*    number of BITMAP_t words in a facet bitmap
*
* double *StoreFacets
*    where facets are stored as block of FacetSize doubles. There is
*    space to store AheadFacets such blocks.
* double *StoreFacetDists
*    for each facet a single double value; MaxFacets blocks.
* BITMAP_t *VertexAdjList, *FacetLiving, *FacetFinal
*    bitmaps with FacetBitmapBlockSize BITMAP_t blocksize. VertexAdjList
*    stores the adjacency list of vertices; there are MaxVertices blocks.
*    Others are facet attributes and they occupy a single block. Each block
*    can store up to MaxFacets bits. FacetLiving tells whether a facet is 
*    among the facets of the current approximation; FacetFinal tells
*     whether the facet belongs to the final polytope.
* BITMAP_t *FacetPositive, *FacetNegative, *FacetWork
*    auxiliary facet bitmaps used during the DD algorithm.
*
* double *FACET_coords(fno)
*    the facet equation; each equation occupies FacetSize doubles, and
*    fno < AheadFacets
* BITMAP_t *FACET_adj(fno)
*    the vertex adjacency list of the facet fno < MaxFacets
* double FACET_dist(fno)
*    the distance associated with facet fno < MaxFacets
*/

static int
  VertexSize,		// size of the vertex coordinate block
  NextVertex,		// next free slot for a vertex
  MaxVertices,		// upper bound for vertex numbers
  VertexBitmapBlockSize,// size of the vertex bitmap block
  FacetSize,		// size of the facet equation block
  NextFacet,		// next free slot for a facet
  MaxFacets,		// number of bits in a facet bitmap
  AheadFacets,		// number of slots when storing facets
  FacetBitmapBlockSize;	// size of the facet bitmap block

/* vertices */
#define StoreVertices	\
    ((double*)memory_slots[M_vertexcoord].ptr)
#define FacetAdjList	\
    ((BITMAP_t*)memory_slots[M_facetadj].ptr)
#define VertexWork	\
    ((BITMAP_t*)memory_slots[M_vertexwork].ptr)
#define VertexArray	\
    ((double*)memory_slots[M_vertexarray].ptr)
#define VERTEX_coord(vno)\
    (StoreVertices+((vno)*VertexSize))
#define VERTEX_adj(vno)	\
    (VertexAdjList+((vno)*FacetBitmapBlockSize))
/* facets */
#define StoreFacets	\
    ((double*)memory_slots[M_facetcoord].ptr)
#define StoreFacetDists	\
    ((double*)memory_slots[M_facetdist].ptr)
#define VertexAdjList	\
    ((BITMAP_t*)memory_slots[M_vertexadj].ptr)
#define FacetLiving	\
    ((BITMAP_t*)memory_slots[M_facetliving].ptr)
#define FacetFinal	\
    ((BITMAP_t*)memory_slots[M_facetfinal].ptr)
#define FacetPositive	\
    ((BITMAP_t*)memory_slots[M_facetpos].ptr)
#define FacetNegative	\
    ((BITMAP_t*)memory_slots[M_facetneg].ptr)
#define FacetWork	\
    ((BITMAP_t*)memory_slots[M_facetwork].ptr)
#define FACET_coords(fno)\
    (StoreFacets+((fno)*FacetSize))
#define FACET_adj(fno)	\
    (FacetAdjList+((fno)*VertexBitmapBlockSize))
#define FACET_dist(fno)	\
    StoreFacetDists[fno]

/***********************************************************************
* Some facet manipulating routines
*
* void move_lastfacet_to(fno)
*   copy the facet at NextFacet to the given place, and decrease
*   NextFacet
*/
#define move_lastfacet_to(fno)	\
  { NextFacet--; \
    memcpy(FACET_coords(fno),FACET_coords(NextFacet),FacetSize*sizeof(double)); \
    memcpy(FACET_adj(fno),FACET_adj(NextFacet),VertexBitmapBlockSize*sizeof(BITMAP_t)); }

/***********************************************************************
* Some auxiliary bitmap functions
*
* is_livingFacet(fno)
*     whether facet fno is an actual facet
* is_finalFacet(fno)
*     whether facet fno is a facet of the final polytope
* clear_Fwork()
*    clear all bits in FacetWork
* copy_Fliving_to_Fwork()
*    copy the bitmap FacetLiving to FacetWork
* clear_facet_in_Fwork(fno)
*    clear bit fno in the bitmap FacetWork
* set_facet_in_Fwork(fno)
*    set bit fno in the bitmap FacetWork
* set_facet_in_Fliving(fno)
*    set bit fno in bitmap FacetLiving
* set_facet_in_Ffinal(fno)
*    set bit fno in bitmap Facetfinal
* intersect_Fwork_with_vertex(vertexno)
*    intersect FacetWork with the adjacency list of the given vertex
* subtract_Fneg_from_Fliving()
*    form the Fliving = Fliving - Fneg bitmap
* clear_PosNegFacets()
*    clear FacetPositive and FacetNegative bitmaps
* add_PositiveFacet(facetno)
*    set bit facetno in FacetPositive
* add_NegativeFacet(facetno)
*    set bot facetno in FacetNegative
*
* copy_VertexWork_to_facet(facetno)
*    copy VertexWork as the adjacency matrix to facet facetno
* clear_FacetAdj_all(facetno)
*    clear the adjacency list of this facet
* set_FacetAdj(facetno,vertexno)
*    set bit vertexno in the adjacency list of facet facetno
* clear_VertexAdj_all(vertexno)
*    clear the adjacency list of a vertex
* set_VertexAdj(vertexno,facetno)
*    set bit facetno in the adjacency list of vertex vertexno
* intersect_VertexAdj_Fliving(vertexno)
*    clear bits in the adjacency list of vertexno which are not in
*    FacetLiving
*/
#define is_livingFacet(idx)	\
	((FacetLiving[(idx)>>packshift]>>((idx)&packmask))&1)
#define is_finalFacet(idx)	\
	((FacetFinal[(idx)>>packshift]>>((idx)&packmask))&1)
#define clear_Fwork()		\
    memset(FacetWork,0,FacetBitmapBlockSize*sizeof(BITMAP_t))
#define copy_Fliving_to_Fwork()	  \
    memcpy(FacetWork,FacetLiving,FacetBitmapBlockSize*sizeof(BITMAP_t))
#define clear_facet_in_Fwork(fno) \
    FacetWork[(fno)>>packshift] &= ~(BITMAP1<<((fno)&packmask))
#define set_facet_in_Fwork(fno)	  \
    FacetWork[(fno)>>packshift] |= BITMAP1<<((fno)&packmask)
#define set_facet_in_Fliving(fno) \
    FacetLiving[(fno)>>packshift] |= BITMAP1<<((fno)&packmask)
#define set_facet_in_Ffinal(fno)  \
    FacetFinal[(fno)>>packshift] |= BITMAP1<<((fno)&packmask)
inline static void intersect_Fwork_with_vertex(int vno)
{int i;
    for(i=0;i<FacetBitmapBlockSize;i++)
        FacetWork[i] &= VERTEX_adj(vno)[i];
}
inline static void subtract_Fneg_from_Fliving(void)
{int i;
    for(i=0;i<FacetBitmapBlockSize;i++)
        FacetLiving[i] &= ~FacetNegative[i];
}

#define clear_PosNegFacets()	  \
    memset(FacetPositive,0,FacetBitmapBlockSize*sizeof(BITMAP_t)); \
    memset(FacetNegative,0,FacetBitmapBlockSize*sizeof(BITMAP_t))
#define add_PositiveFacet(fno)	  \
    FacetPositive[(fno)>>packshift] |= BITMAP1<<((fno)&packmask)
#define add_NegativeFacet(fno)	  \
    FacetNegative[(fno)>>packshift] |= BITMAP1<<((fno)&packmask)
#define copy_VertexWork_to_facet(fno) \
    memcpy(FACET_adj(fno),VertexWork,VertexBitmapBlockSize*sizeof(BITMAP_t))
#define clear_FacetAdj_all(fno)	  \
    memset(FACET_adj(fno),0,VertexBitmapBlockSize*sizeof(BITMAP_t))
#define set_FacetAdj(fno,vno)	  \
    FACET_adj(fno)[(vno)>>packshift] |= BITMAP1 << ((vno)&packmask)
#define clear_VertexAdj_all(vno)  \
    memset(VERTEX_adj(vno),0,FacetBitmapBlockSize*sizeof(BITMAP_t))
#define set_VertexAdj(vno,fno)	  \
    VERTEX_adj(vno)[(fno)>>packshift] |= BITMAP1 << ((fno)&packmask)
inline static void intersect_VertexAdj_Fliving(int vno)
{int i;
    for(i=0;i<FacetBitmapBlockSize;i++){
        VERTEX_adj(vno)[i] &= FacetLiving[i];
    }
}
/***********************************************************************
* Statistics
* The struct dd_stats contains statistics on the algorithm.
*
* get_dd_facetno()
*   computes the number of living and final facets.
*/
DD_STATS dd_stats;

void get_dd_facetno(void) // get the number living and final facets
{int fno,i; BITMAP_t v;
    dd_stats.living_facets_no=-1;
    dd_stats.final_facets_no=-1;
    fno=0; for(i=0;i<FacetBitmapBlockSize;i++){
        v=FacetLiving[i];
        if(v==~BITMAP0){ dd_stats.living_facets_no+=(1<<packshift); }
        else while(v){
            if(v&1) dd_stats.living_facets_no++;
            v>>=1;
        }
        v=FacetFinal[i];
        if(v==~BITMAP0){ dd_stats.final_facets_no +=(1<<packshift); }
        else while(v){
            if(v&1) dd_stats.final_facets_no++;
            v>>=1;
        }
        fno += (1<<packshift);
    }
}
/***********************************************************************
* Interrogate facets and vertices
*
* int get_next_facet(from)
*   return the next living but not final facet number starting at
*   from. If no such vertices are, return -1.
*   When from=-1 and RandomFacet is set, pick the starting number
*   randomly.
*
* void mark_facet_as_final(fno)
*   mark the facets as final, i.e. facet of the final polytope.
*
* int vertex_num(void)
*   return the number of vertices added so far
*
* int facet_num(void)
*   return the number of the facets of the actual approximating polytope;
*   this is the number of facets of the final polytope when the algorithm
*   terminates
*
* void get_facet_into(fno, double v[0:dim])
*   store the facet equation to the provided space. The equation is scaled
*   and rounded to the closet integer when the difference is small
*
* void print_vertex(report_type, double dir, double coords[0:DIM-1](
*   report the vertex using fractional or floating format.
*
* void print_vertices(report_type where)
*   report all stored vertices using print_vertex().
*
* void print_facets(report_type where)
*   report all facets using scaling and rounding format.
*
* int mrandom(v)
*   return an integer between 0 and v-1, approximately uniformly.
*/

/** random integer between 0 and v-1 **/
static inline int mrandom(int v)
{ return v<=1?0 : (random()&0x3fffffff)%v; }

/** return the next living facet from 'from' or from a random location **/
int get_next_facet(int from)
{int fno,i,j; BITMAP_t v;
    dd_stats.facetenquires++;
    if(from<0){
        if(PARAMS(RandomFacet)){ // start from a random place
            dd_stats.facetenquires--;
            fno=get_next_facet(mrandom(NextFacet));
            if(fno>=0) return fno;
        }
        from=0; 
    }
    else if((from>>packshift) >= FacetBitmapBlockSize)
    { return -1; }
    j=from & packmask;
    fno = from-j;
    i = fno>>packshift;
    if(j){
        v=(FacetLiving[i] & ~FacetFinal[i])>>j;
        while(v){
            if(v&1){ return fno+j; }
            j++; v>>=1;
        }
        fno += (1<<packshift); i++;
    }
    for(j=0;i<FacetBitmapBlockSize;i++){
        v= FacetLiving[i] & ~FacetFinal[i];
        while(v){
            if(v&1){ return fno+j; }
            j++; v>>=1;
        }
        fno += (1<<packshift);
    }
    return -1; 
}

/** mark the facet as final **/
void mark_facet_as_final(int fno)
{    FacetFinal[fno>>packshift] |= BITMAP1 << (fno&packmask); }

/** number of added vertices (except the ideal ones) **/
int vertex_num(void)
{ return NextVertex-DIM; }

/** number of living, non-ideal facets **/
int facet_num(void)
{int i,j;
    j=0; for(i=1;i<NextFacet;i++)if(is_livingFacet(i)) j++;
    return j;
}

/** round x to an integer if closer than 3e-9 **/
#define SCALE_EPS	PARAMS(ScaleEps)
static double round(double x)
{double res,y;
     if(x<0.0){
         res=(int)(-x+0.5); y=x+res;
         if(y<SCALE_EPS && y>-SCALE_EPS) return -res;
     } else {
         res=(int)(x+0.5); y=x-res;
         if(y<SCALE_EPS && y>-SCALE_EPS) return res;
     }
     return x;
}

/** check if the argument is close to a (small) integer value **/
inline static int closetoint(double x)
{
    if(x<0.0) x=-x;
    x -= (int)(x+0.5); 
    return (x<1e-8 && x>-1e-8);
}

/** gcd, lcm: compute gcd and lcm of two integers **/
static int gcd(int a, int b)
{   if(a==0 || a==b) return b;
    if(b==0) return a;
    if(a<b) return gcd(b%a,a);
    return gcd(a%b,b);
}
static int lcm(int a,int b)
{   if(b==1 || a==b) return a;
    if(a==1) return b;
    return a*(b/gcd(a,b));
}

/** check if x is close to a rational with denominator at most 1000 */
static int denum(double x)
{int i;
    for(i=1;i<1000;i++)if(closetoint(i*x)) return i;
    return 1;
}

/** store the scaled and rounded facet equation at the given location **/
void get_facet_into(int fno, double *v)
{int j,d;
    d=1; for(j=0;d<130000 && j<=DIM;j++)d=lcm(d,denum(FACET_coords(fno)[j]));
    for(j=0;j<=DIM;j++)v[j]=round(d*FACET_coords(fno)[j]);
}

/** formatted output of a vertex coordinate **/
static char *formatvalue(double v)
{static char buf[80]; int d;
    if(!closetoint(v) && PARAMS(VertexAsFraction)){
        d=denum(v); if(closetoint(d*v)){
            sprintf(buf,"%d/%d",(int)(round(d*v)),d);
            return buf;
        }
    }
    sprintf(buf,"%.14g",round(v));
    return buf;
}

/** print a single vertex in fractional or float format **/
void print_vertex(report_type where, double dir, const double v[/*0:DIM-1*/])
{int j;
    for(j=0;j<PARAMS(ProblemObjects);j++){ // DIM might not be set
        report(where," %s",formatvalue(dir*v[j]));
    }
    report(where,"\n");
}

/** print all vertices **/
void print_vertices(report_type where)
{int i; double dir;
    dir = PARAMS(Direction) ? -1.0 : +1.0;
    for(i=DIM;i<NextVertex;i++){
        report(where,"V ");
        print_vertex(where,dir,VERTEX_coord(i));
    }
}

/** print scaled and rounded facet equations **/
void print_facets(report_type where)
{int i,j,d;
//    report(where,"C Total number of facets: %d\n",facet_num());
    for(i=1;i<NextFacet;i++)if(is_livingFacet(i)){
        report(where,"F ");
        d=1; for(j=0;d<130000 &&j<=DIM;j++)d=lcm(d,denum(FACET_coords(i)[j]));
        for(j=0;j<=DIM;j++)
           report(where," %.14lg",round(d*FACET_coords(i)[j]));
        report(where,"\n");
    }
}

/***********************************************************************
* initialize all data structures and add the very first vertex
*
* int VERTEXARRAY_STEPSIZE
*   the step size of VertexArray which stores vertices.
*
* int init_dd(int dim, double coords[0:dim-1])
*   dim:    the dimension of the problem; should be not too large
*   coords: coordinates of the very first vertex.
* return value
*   0:  OK
*   1:  some error: dimension is too large, or not enough memory
*/

#define VERTEXARRAY_STEPSIZE	(3*DIM)

int init_dd(int dimension, double *coords)
{int i,j;
    DIM=dimension; // set the dimension
    if(DIM<1||DIM>MAXIMAL_ALLOWED_DIMENSION){
        report(R_fatal,"init_dd: dimension %d is out of allowed range (1..%d)\n",
           dimension,MAXIMAL_ALLOWED_DIMENSION);
           return 1;
    }
    MaxVertices = DD_INITIAL_VERTEXNO;
    AheadFacets = MaxFacets = DD_INITIAL_FACETNO;
    while(MaxVertices+5<DIM) MaxVertices += (DD_VERTEX_ADDBLOCK<<packshift);
    // how many doubles store a facet and a vertex
    FacetSize=DIM+1; VertexSize=DIM;
    // bitmaps, allocate managed memory
    VertexBitmapBlockSize = (MaxVertices+packmask)>>packshift;
    FacetBitmapBlockSize  = (MaxFacets+packmask)>>packshift;
    // clear statistics data
    memset(&dd_stats,0,sizeof(dd_stats));
    dd_stats.vertices_allocated_no=1;
    dd_stats.vertices_allocated=MaxVertices;
    dd_stats.facets_allocated_no=1;
    dd_stats.facets_allocated=MaxFacets;
    yalloc(BITMAP_t,M_facetadj,AheadFacets,VertexBitmapBlockSize); // FacetAdjList
    yalloc(BITMAP_t,M_vertexwork,1,VertexBitmapBlockSize); // VertexWork
    yalloc(BITMAP_t,M_facetfinal,1,FacetBitmapBlockSize); // FacetFinal
    yalloc(BITMAP_t,M_facetliving,1,FacetBitmapBlockSize); // FacetLiving
    yalloc(BITMAP_t,M_facetpos,1,FacetBitmapBlockSize);  // FacetPositive
    yalloc(BITMAP_t,M_facetneg,1,FacetBitmapBlockSize); // FacetNegative
    yalloc(BITMAP_t,M_facetwork,1,FacetBitmapBlockSize); // FacetWork
    yalloc(BITMAP_t,M_vertexadj,MaxVertices,FacetBitmapBlockSize); // VertexAdjList
    yalloc(double,M_vertexcoord,MaxVertices,VertexSize); // StoreVertices
    yalloc(double,M_facetdist,AheadFacets,1); // StoreFacetDists
    yalloc(double,M_facetcoord,AheadFacets,FacetSize); // StoreFacets
    yalloc(double,M_vertexarray,VERTEXARRAY_STEPSIZE,DIM+1); // VertexArray
    if(OUT_OF_MEMORY) return 1;
    dd_stats.memory_allocated_no=1;
    // vertices and facets
    NextVertex=0;
    NextFacet=0;
    /** create the ideal vertices **/
    for(i=0;i<DIM;i++){ // the i-th ideal vertex
        VERTEX_coord(NextVertex)[i]=1.0;
        clear_VertexAdj_all(NextVertex);
        set_VertexAdj(NextVertex,0); // it is on the ideal facet
        // and on all facets except the i-th one
        for(j=0;j<DIM;j++) if(i!=j)set_VertexAdj(NextVertex,j+1);
        NextVertex++;
    }
    // the very first real vertex
    for(j=0;j<DIM;j++) VERTEX_coord(NextVertex)[j]=coords[j];
    clear_VertexAdj_all(NextVertex);
    for(j=0;j<DIM;j++)set_VertexAdj(NextVertex,j+1);
    NextVertex++;
    // the ideal facet
    FACET_coords(NextFacet)[DIM]=1.0;
    clear_FacetAdj_all(NextFacet);
    for(j=0;j<DIM;j++)set_FacetAdj(NextFacet,j);
    set_facet_in_Fliving(NextFacet); set_facet_in_Ffinal(NextFacet);
    NextFacet++;
    // other facets
    for(i=0;i<DIM;i++){
        FACET_coords(NextFacet)[i]=1.0;
        FACET_coords(NextFacet)[DIM]=-coords[i];
        clear_FacetAdj_all(NextFacet);
        set_FacetAdj(NextFacet,DIM);
        for(j=0;j<DIM;j++) if(i!=j) set_FacetAdj(NextFacet,j);
        set_facet_in_Fliving(NextFacet);
        NextFacet++;
    }
    return 0;
}

/***********************************************************************
* Request memory to accommodate more vertices or more facets.
*
* void allocate_vertex_block()
*   add space for DD_VERTEX_ADDBLOCK<<packshift more vertices, then
*    expand StoreVertices, VertexAdjList;
*    adjust VertexBitmapBlockSize by DD_VERTEX_ADDBLOCK
*    reallocate all vertex bitmaps (FacetAdjList and WertexWork)
*    Called at the very beginning of a new iteration when necessary.
*
* void allocate_facet_block()
*   add space for DD_FACET_ADDBLOCK<<packshift more facets. Only
*     increase storage place in StoreFacets and FacetAdjList, but
*     do not change facet bitmaps.
*
* void expand_facet_bitmaps()
*   expand all facet bitmaps to a size which can hold AheadFacets
*     many bits. Adjust MaxFacets accordingly. Keep the content of
*     facet bitmaps Final, Living, Work; other can be dropped.
*/

static void allocate_vertex_block(void)
{   // extend Vertices and FacetBitmap to accommodate more stuff;
    dd_stats.vertices_allocated_no ++;
    dd_stats.vertices_allocated += (DD_VERTEX_ADDBLOCK<<packshift);
    MaxVertices += (DD_VERTEX_ADDBLOCK<<packshift);
    VertexBitmapBlockSize += DD_VERTEX_ADDBLOCK;
    // tell the memory handling part how much storage space we would need.
    yrequest(double,M_vertexcoord,MaxVertices,VertexSize);
    yrequest(BITMAP_t,M_vertexadj,MaxVertices,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_vertexwork,1,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_facetadj,AheadFacets,VertexBitmapBlockSize);
    if(reallocmem()){  // out of memory, don't increase the values
        MaxVertices -= (DD_VERTEX_ADDBLOCK<<packshift);
        VertexBitmapBlockSize -= DD_VERTEX_ADDBLOCK;
    }
}

static void allocate_facet_block(void)
{   if(OUT_OF_MEMORY) return;
    dd_stats.facets_allocated_no ++;
    dd_stats.facets_allocated += DD_FACET_ADDBLOCK<<packshift;
    AheadFacets += DD_FACET_ADDBLOCK<<packshift;
    yrequest(double,M_facetcoord,AheadFacets,FacetSize);
    yrequest(BITMAP_t,M_facetadj,AheadFacets,VertexBitmapBlockSize);
    if(reallocmem()){ // out of memory
        AheadFacets -= DD_FACET_ADDBLOCK<<packshift;
    }
}

static void expand_facet_bitmaps(void)
{int newblocksize;
    dd_stats.expand_facetbitmaps_no ++;
    newblocksize=(AheadFacets+packmask)>>packshift;
    yrequest(double,M_facetdist,AheadFacets,1);
    yrequest(BITMAP_t,M_facetfinal,1,newblocksize);
    yrequest(BITMAP_t,M_facetliving,1,newblocksize);
    yrequest(BITMAP_t,M_facetpos,1,newblocksize);
    yrequest(BITMAP_t,M_facetneg,1,newblocksize);
    yrequest(BITMAP_t,M_facetwork,1,newblocksize);
    yrequest(BITMAP_t,M_vertexadj,MaxVertices,newblocksize);
    if(reallocmem()){ // out of memory
        return;
    }
    FacetBitmapBlockSize = newblocksize;
    MaxFacets=AheadFacets;
}

/***********************************************************************
* Get the next facet number
* int get_new_facetno()
*    return the new facet number, asking for memory if necessary.
*    The return value is -1 if cannot allocate more memory.
*/
inline static int get_new_facetno(void)
{int i;
    i=NextFacet; 
    if(NextFacet>=AheadFacets){
      allocate_facet_block();
      if(OUT_OF_MEMORY) return -1;
    }
    NextFacet++; return i;
}

/***********************************************************************
* Compute facet equation from the vertices it is adjacent to.
*
* int solve_lineq(int d)
*   solve the system of homogeneous linear equation stored in
*   VertexArray[dim+1,d]. Use Gauss elimination: get the largest value
*   in the first column A[0,], swap this row and the first row; and
*   subtract from all rows so that the first column will be all zero,
*   etc. The matrix rank should be exactly dim (so that we have a
*   non-trivial solution). The solution is returned in A[0..dim].
*   PARAMS(LineqEps) is the rank threshold.
*   Returns 1 if the matrix rank is not dim; otherwise returns zero.
*
* void recalculate_facet_eq(int fno)
*   calculate facet equation form the vertices adjacent to it. Complain
*   if out of memory, the system is degenerate, or the old and new coeffs
*   are too far from each other.
*
* void recalculate_facets(void)
*   go over all facets, and recalculate their equations.
*/
static int solve_lineq(int d)
{int D1,i,j,jmax,col,zerocol; double v,vmax;
    D1=DIM+1; zerocol=-1;
#define A(i,j)	VertexArray[(i)*D1+(j)]
    for(col=0,j=0;col<=DIM;j++,col++){
       /* get the largest value of column col to A[j,col] */
       jmax=j; vmax=0.0;
       for(i=j;i<d;i++){
           v=A(i,col); if(v<0.0) v=-v;
           if(vmax<v){jmax=i;vmax=v;}
       }
       if(vmax< PARAMS(LineqEps)){ /* too small value */
           if(zerocol>=0){
              report(R_err,"solve_lineq: rank is too small (zerocol=%d, col=%d, max=%lg)\n"
                 ,zerocol,col,vmax);
              return 1; /* error */
           }
           zerocol=col; j--;
           continue;
       }
       if(jmax!=j){ /* swap rows jmax and j */
          for(i=col;i<=DIM;i++){
             v=A(j,i); A(j,i)=A(jmax,i); A(jmax,i)=v;
          }
       }
       v=1.0/A(j,col);
       for(i=col+1;i<=DIM;i++){ A(j,i)*=v; }
       A(j,col)=1.0;
       /* and subtract row j from all rows */
       for(jmax=0;jmax<d;jmax++) if(jmax!=j){
           v=A(jmax,col);
           for(i=col+1;i<=DIM;i++){ A(jmax,i) -= A(j,i)*v; }
           A(jmax,col)=0.0;
       }
    }
    if(zerocol<0){
        report(R_err,"solve_lineq: rank is too large, tolerance LineqEps=%lg\n",
                     PARAMS(LineqEps));
        return 1; /* error */
    }
    /* compute values to row 0 */
    for(j=0,col=0;col<=DIM;j++,col++){
        if(col==zerocol){j--; continue; }
        A(0,col) = round(-A(j,zerocol));
    }
    A(0,zerocol)=1.0;
    return 0;
#undef A
}

static void recalculate_facet_eq(int fno)
{double v,s; int i,j,vno; int amax,an; BITMAP_t fc;
#define A(i,j)	VertexArray[(i)*(DIM+1)+(j)]
    // collect all vertices adjacent to facet fno
    amax=VERTEXARRAY_STEPSIZE; /* we have at least that many slots */
    an=0; vno=0;
    for(i=0;i<VertexBitmapBlockSize;i++){
        j=vno; fc=FACET_adj(fno)[i];
            while(fc){
            if(fc&1){ /* store vertex j to A */
                if(an>=amax){
                    amax+=VERTEXARRAY_STEPSIZE;
                    yrequest(double,M_vertexarray,amax,DIM+1);
                    if(reallocmem()) return; // out of memory
                }
                memcpy(&(A(an,0)),VERTEX_coord(j),DIM*sizeof(double));
                A(an,DIM)= j<DIM ? 0.0 : 1.0 ;
                an++;
            }
            j++; fc>>=1;
        }
        vno += (1<<packshift);
    }
    // all vertices are collected, compute facet equation
    if(an<DIM){
        report(R_err,"recalculate_facet: facet %d has only %d(<DIM=%d) adjacent vertices\n",fno,vno,DIM);
        dd_stats.numerical_error++;
        return;
    }
    if(solve_lineq(an)){
        report(R_err,"recalculate_facet: adjacency vertex list of facet %d is degenerate\n",fno);
        dd_stats.numerical_error++;
        return;
    }
    s=0.0; for(i=0;i<DIM;i++) s+= VertexArray[i]; s=1.0/s;
    for(i=0;i<=DIM;i++){
        VertexArray[i]*=s;
        v=FACET_coords(fno)[i]-VertexArray[i];
        if(v>PARAMS(FacetRecalcEps) || v < -PARAMS(FacetRecalcEps)){
            report(R_warn,"recalculate_facet: numerical instability at facet %d, coord %d (%lg)\n", fno,i,v);
            dd_stats.instability_warning++;
        }
        FACET_coords(fno)[i]=VertexArray[i];
    }
#undef A
}

void recalculate_facets(void)
{int fno; 
    for(fno=1;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        recalculate_facet_eq(fno);
    }
}

/***********************************************************************
* double vertex_distance(double *coords, int fno)
*   compute the distance of the vertex given as the first argument
*   from the facet given with its number
*
* int store_vertex(double *coords)
*   check if the vertex is new, and if yes, store it, and disable
*   further vertex processing. Return 1 if this vertex was not seen
*   before; and zero otherwise.
*
* double probe_vertex(double *coords)
*   return the number of ridge tests for this vertex.
*/

#define DD_EPS_EQ	PARAMS(PolytopeEps)

inline static double vertex_distance(double *coords, int fno)
{double d=0.0; int i; double *fcoords;
    fcoords=FACET_coords(fno);
    for(i=0;i<DIM;i++){
       d+= (*coords)*(*fcoords);
       coords++; fcoords++;
    }
    d += (*fcoords);
    FACET_dist(fno)=d;
    return d;
}

int store_vertex(double *coords) /* store vertex if new */
{int i,j; int thisvertex; double d;
    for(i=DIM;i<NextVertex;i++){ /* check if this is a new vertex */
        for(j=0;j<DIM;j++){
            d=VERTEX_coord(i)[j]-coords[j];
            if(d>DD_EPS_EQ || d < -DD_EPS_EQ) break;
        }
        if(j==DIM) return 0; /* this is the same as the one stored */
    }
    dd_stats.iterations++;
    if(NextVertex>=MaxVertices){
        allocate_vertex_block();
        if(OUT_OF_MEMORY) return 1; // out of memory
    }
    thisvertex=NextVertex; NextVertex++;
    for(i=0;i<DIM;i++){ VERTEX_coord(thisvertex)[i]=coords[i]; }
    clear_VertexAdj_all(thisvertex); // clear adjacency list
    return 1; /* new vertex */
}

double probe_vertex(double *coords) /* return the number of ridge tests */
{double d; int fno,posfacets,negfacets;
    dd_stats.probevertex++;
    posfacets=0; negfacets=0;
    for(fno=0;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        d=vertex_distance(coords,fno);
        if(d>DD_EPS_EQ){ // positive side
            posfacets++;
        } else if(d<-DD_EPS_EQ){ // negative
            negfacets++;
        }
    }
    return (double)posfacets * (double)negfacets;
}

/***********************************************************************
* Routines for vertex enumeration
*
* int facet_intersection(f1,f2)
*    intersect the vertex adjacency lists of f1 and f2, store it in
*    VertexWork. Return the number of vertices adjacent to both
*    f1 and f2.
*
* int is_ridge(f1,f2)
*    combinatorial test to check whether the intersection of f1 and
*    f2 forms a ridge. VertexWork contains vertices adjacent to both
*    f1 and f2. If there are <DIM-1 such vertices, it is not a ridge.
*    Otherwise compute the intersection of living facets and the facet
*    adjacency lists of vertices adjacent to both f1 and f2. This is
*    a ridge if and only if the intersection contains f1 and f2 only.
*
* void create_new_facet(f1,f2,vno)
*    create a new facet through the ridge formed by f1 and f2, and
*    the given vertex. f1*v is negative, f2*v is positive, and these
*    values are stored in FACET_dist(). Recalculate the equation when
*    ExactFacetEq parameter is set
*
* void search_ridges_with(f1,vno)
*    Vertex vno is on the negative side of f1; go over all positive
*    facets and call create_new_facet()
*
* void make_facet_living(fno)
*    set the "living" flag for this facet, and add it to the adjacency
*    list of all vertices it is adjacent to.
*
* void add_new_vertex(double vertex[0:dim-1])
*    add a new vertex which is outside the convex hull of the present
*    approximation. Split facets into positive, negative and zero sets
*    depending where the new vertex is. For each pair of positive/
*    negative facets, check if it is a ridge; If yes, add a new facet.
*    The throw away negative facets, and add the newly created facets
*    to the approximation.
*/

/** compute how many vertices are adjacent to both f1 and f2 **/
inline static int facet_intersection(int f1, int f2)
{int i,total; register BITMAP_t *L1,*L2;
#ifdef BITMAP_32
    BITMAP_t v;
    total=0; L1=FACET_adj(f1); L2=FACET_adj(f2);
    for(i=0;i<VertexBitmapBlockSize;i++,L1++,L2++){
        v=VertexWork[i]= (*L1)&(*L2);
        while(v){if(v&1) total++; v>>=1; }
    }
#else
    SPLIT_BITMAP_t vv;
    total=0; L1=FACET_adj(f1); L2=FACET_adj(f2);
    for(i=0;i<VertexBitmapBlockSize;i++,L1++,L2++){
        vv.u=VertexWork[i]= (*L1)&(*L2);
        while(vv.i.i1){if(vv.i.i1&1) total++; vv.i.i1>>=1; }
        while(vv.i.i2){if(vv.i.i2&1) total++; vv.i.i2>>=1; }
    }
#endif
    return total;
}

/** combinatorial test whether f1 and f2 intersect in a ridge **/
inline static int is_ridge(int f1,int f2)
{int vertexno,i,j;
     // get vertices adjacent to both f1 and f2 into VertexWork
    if(facet_intersection(f1,f2) < DIM-1){
         return 0; // no
    }
    copy_Fliving_to_Fwork();
#ifdef BITMAP_32
    {BITMAP_t v;
    vertexno=0;
    for(i=0;i<VertexBitmapBlockSize;i++){
        j=vertexno; v=VertexWork[i];
        while(v){
            if(v&1){ intersect_Fwork_with_vertex(j); }
            j++; v>>=1;
        }
        vertexno += (1<<packshift);
    }}
#else
    {uint32_t *vv, v;
    vertexno=0; vv=(uint32_t *)VertexWork;
    for(i=0;i<VertexBitmapBlockSize;i++){
        j=vertexno; v=*vv;
        while(v){
            if(v&1){ intersect_Fwork_with_vertex(j); }
            j++; v>>=1;
        }
        vertexno += (1<<(packshift-1)); j=vertexno; vv++; v=*vv; 
        while(v){
            if(v&1){ intersect_Fwork_with_vertex(j); }
            j++; v>>=1;
        }
        vv++; vertexno += (1<<(packshift-1));
    }}
#endif
    // only f1 and f2 should remain, if any other remains, no
    clear_facet_in_Fwork(f1); clear_facet_in_Fwork(f2);
    for(i=0;i<FacetBitmapBlockSize;i++){
        if(FacetWork[i]){ 
            return 0; // no
        }
    }
    return 1; // yes
}

/** create a new facet through the ridge f1, f2 and vertex vno **/
inline static void create_new_facet(int f1, int f2, int vno)
{int newf; double d1,d2,d; int i;
    newf=get_new_facetno();
    if(newf<0) return; // no memory
    copy_VertexWork_to_facet(newf); set_FacetAdj(newf,vno);
    // compute the coefficients, f1<0, f2 >0
    if(f2==0){ // ideal facet
        for(i=0;i<=DIM;i++){
            FACET_coords(newf)[i]=FACET_coords(f1)[i];
        }
        FACET_coords(newf)[DIM] -= FACET_dist(f1);
    } else {   // f1(v)=-d1, f2(v)=d2; (d2*f1 + d1*f2)(v)=0.0
        d1=-FACET_dist(f1); d2=FACET_dist(f2);
        d=1.0/(d1+d2); d1 *= d; d2 *= d;
        for(i=0;i<=DIM;i++){
            FACET_coords(newf)[i]=d2*FACET_coords(f1)[i]+d1*FACET_coords(f2)[i];
        }
    }
    if(PARAMS(ExactFacetEq)) recalculate_facet_eq(newf);
}

/** facet f1 is negative; go over all positive facets **/
inline static void search_ridges_with(int f1, int newvertex)
{int i,j,f2; BITMAP_t f2c;
    f2=0;
    for(i=0;i<FacetBitmapBlockSize;i++){
        if(OUT_OF_MEMORY) return;
        j=f2; f2c=FacetPositive[i];
        while(f2c){
            if((f2c&1) && is_ridge(f1,j)){
                create_new_facet(f1,j,newvertex);
            }
            j++; f2c>>=1;
        }
        f2+=(1<<packshift);
    }
}

/** make the facet living, adjust the adjacency list of its vertices **/
static void make_facet_living(int fno)
{int vno,i,j; BITMAP_t fc;
    set_facet_in_Fliving(fno);
    vno=0; for(i=0;i<VertexBitmapBlockSize;i++){
        j=vno; fc=FACET_adj(fno)[i];
        while(fc){
            if(fc&1){set_VertexAdj(j,fno);}
            j++; fc>>=1;
        }
        vno+=(1<<packshift);
    }
}

/** add a new vertex to the approximation **/
void add_new_vertex(double *coords)
{int i,j,fno,vno; int thisvertex; BITMAP_t fc; double d;
 int posfacets,negfacets,zerofacets,newfacet;
    dd_stats.iterations++;
    if(NextVertex >= MaxVertices){
        allocate_vertex_block();
        if(OUT_OF_MEMORY) return;
    }
    thisvertex=NextVertex; NextVertex++;
    for(i=0;i<DIM;i++){ VERTEX_coord(thisvertex)[i]=coords[i]; }
    clear_VertexAdj_all(thisvertex); // clear the adjacency list
    clear_PosNegFacets();
    // split facets into positive, negative and zero parts
    posfacets=0; negfacets=0; zerofacets=0;
    for(fno=0;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        d=vertex_distance(coords,fno);
        if(d>DD_EPS_EQ){ // positive side 
            add_PositiveFacet(fno);
            posfacets++;
        } else if(d<-DD_EPS_EQ){ // negative
            add_NegativeFacet(fno);
            negfacets++;
        } else { // this is adjacent to our new vertex
            set_VertexAdj(thisvertex,fno);
            set_FacetAdj(fno,thisvertex);
            zerofacets++;
        }
    }
    // add statistics: posfacets*negfacets
    d=(double)posfacets*(double)negfacets;
    if(dd_stats.max_tests<d)dd_stats.max_tests=d;
    dd_stats.avg_tests = (dd_stats.iterations-1)*dd_stats.avg_tests*(1.0/dd_stats.iterations)
      + (1.0/dd_stats.iterations)*d;
    // this is the critical part: for each positive/negative facet pair
    // compute whether this is a new facet.
    newfacet=NextFacet; // how many new facets are added at this step
    fno=0;for(i=0;i<FacetBitmapBlockSize;i++){
        j=fno; fc=FacetNegative[i];
        while(fc){
            if(fc&1){search_ridges_with(j,thisvertex);}
            j++; fc>>=1;
        }
        fno+=(1<<packshift);
    }
    if(OUT_OF_MEMORY || dobreak){ NextFacet=newfacet; return; }
    // new facets are from newfacet ... NextFacet.
    // calculate the number of facets of the new polytope
    zerofacets += posfacets + (NextFacet-newfacet);
    if(dd_stats.max_facets < zerofacets) dd_stats.max_facets = zerofacets;
    dd_stats.last_facets[dd_stats.iterations & 15]=zerofacets;
    if(dd_stats.max_facetsadded<NextFacet-newfacet){ dd_stats.max_facetsadded=NextFacet-newfacet; }
    dd_stats.avg_facetsadded = (dd_stats.iterations-1)*dd_stats.avg_facetsadded*(1.0/dd_stats.iterations)
         + (1.0/dd_stats.iterations)*(NextFacet-newfacet);
    // delete negative facets from FacetLiving
    subtract_Fneg_from_Fliving();
    if(NextFacet <= MaxFacets){
        // add new facets to the living list
        while(newfacet<NextFacet){ // distance must be zero ...
            make_facet_living(newfacet);
            newfacet++;
        }
        return;
    }
    dd_stats.facets_compressed_no ++;
    for(vno=0;vno<NextVertex;vno++){ // clear the adjacency list of vertices
        intersect_VertexAdj_Fliving(vno);
    }
    clear_Fwork(); // collect new facets in Fwork
    // move new facets into free facet slots
    // Do it backward starting at NextFacet-1
    fno=0; for(i=0;i<FacetBitmapBlockSize;i++){
        j=fno; fc=~FacetLiving[i]; // complement ...
        while(fc){
            if((fc&1)){
                if(j>=newfacet) goto finish_compress;
                move_lastfacet_to(j); set_facet_in_Fwork(j);
                if(NextFacet<=newfacet) goto finish_compress;
            }
            j++; fc>>=1;
        }
        fno += (1<<packshift);
    }
  finish_compress:
    if(NextFacet>MaxFacets){
        expand_facet_bitmaps();
        // if no memory, throw away some newly generated facets
        if(OUT_OF_MEMORY) NextFacet=MaxFacets;
    }
    // from newfacet until NextFacet indicate that they are new facets
    while(newfacet<NextFacet){
        set_facet_in_Fwork(newfacet); newfacet++;
    }
    // merge FacetWork to FacetLiving
    fno=0; for(i=0;i<FacetBitmapBlockSize;i++){
        j=fno; fc=FacetWork[i];
        while(fc){
            if(fc&1){ make_facet_living(j); }
            j++; fc>>=1;
        }
        fno += (1<<packshift);
    }
    // move back NextFacet as far as possible ...
    while(NextFacet>0 && !is_livingFacet(NextFacet-1)) NextFacet--;
}

/*=====================================================================*/
/** consistency check: vertex bitmaps are set OK **/
int check_consistency(void)
{int i,vno,fno,errno; double w,*f,*v;
    errno=0;
    // check that living facets have >=0 coords and they add up to 1.0 */
    for(fno=1;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        w=-1.0; f=FACET_coords(fno);
        for(i=0;i<DIM;i++,f++){
            w += (*f); if(*f<-DD_EPS_EQ){
                report(R_err,"Consistency error 1: facet %d has negative coeff (%lg)\n",fno,*f);
                errno++;
            }
        }
        if(w>DD_EPS_EQ || w<-DD_EPS_EQ){
            report(R_err,"Consistency error 2: sum of coeffs of facet %d differs from 1.0 (%lg)\n",fno,w);
            errno++;
        }
    }
    // go over all vertices
    for(vno=DIM;vno<NextVertex;vno++){
        for(fno=1;fno<NextFacet;fno++)if(is_livingFacet(fno)){
            // check if vno is adjacent to fno
            f=FACET_coords(fno); v=VERTEX_coord(vno); w=0.0;
            for(i=0;i<DIM;i++,f++,v++){ w += (*f)*(*v); }
            w += (*f); if(w<-DD_EPS_EQ){ 
               report(R_err,"Consistency error 3: vertex %d is on the negative side of facet %d\n",vno,fno);
               errno++;
            }
            if(w>DD_EPS_EQ){ // not adjacent
                if(1&(FACET_adj(fno)[vno>>packshift]>>(vno&packmask))){
                    report(R_err,"Consistency error 4: vertex %d and facet %d are NOT adjacent\n",vno,fno);
                    errno++;
                }
                if(1&(VERTEX_adj(vno)[fno>>packshift]>>(fno&packmask))){
                    report(R_err,"Consistency error 5: vertex %d and facet %d are NOT adjacent\n",vno,fno);
                    errno++;
                }
            } else { // adjacent
                if(!(1&(FACET_adj(fno)[vno>>packshift]>>(vno&packmask)))){
                    report(R_err,"Consistency error 6: vertex %d and facet %d are adjacent\n",vno,fno);
                    errno++;
                }
                if(!(1&(VERTEX_adj(vno)[fno>>packshift]>>(fno&packmask)))){
                    report(R_err,"Consistency error 7: vertex %d and facet %d are adjacent\n",vno,fno);
                    errno++;
                }
            }
        }
    }
    return errno;
}

/* EOF */

