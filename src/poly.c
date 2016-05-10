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
* Content in the main slots are needed for the duration of the algorithm,
* and the blocksize can change. Temporary slots are used within a single
* iteration, and blocksize is fixed.
*
* int M_RESERVEDSLOTS
*    memory slots reserved for threads. Each thread uses temporary slots
*    vertexlist, facetwork, vertexarray, newfacetcoord, newfacetadj,
*    thus this should be 5*(max_threads-1) */
#define M_RESERVEDSLOTS	0 /* number of reserved slots */

typedef enum {	/* main memory slots */
M_vertexcoord	= 0,	/* vertex coordinates; StoreVertices */
M_facetcoord,		/* facet coordinates; StoreFacets */
M_vertexadj,		/* adjacency list of vertices; VertexAdjList*/
M_facetadj,		/* adjacency list of facets; FacetAdjList */
M_facetliving,		/* facet bitmap; FacetLiving */
M_facetfinal,		/* facet bitmap; FacetFinal */
M_MAINSLOTS,		/* last main slot index */
		/* temporary memory slots */
M_facetdist=M_MAINSLOTS,/* facet distances; StoreFacetDists */
M_facetposneg,		/* indices of positive/negative facets; FacetPosnegList */
M_vertexlist,		/* vertex indices in the intersection of two facets */
M_facetwork,		/* facet bitmap; FacetWork */
M_vertexarray,		/* calculating facet equations; VertexArray */
M_newfacetcoord,	/* new facet coordinates */
M_newfacetadj,		/* adjacency list of new facets */
		/* reserved slots used by threads */
M_reserved,
		/* total number of managed memory */
M_MSLOTSTOTAL	= M_reserved + M_RESERVEDSLOTS
} memslot_t;

typedef struct { /* memory slot structure */
  size_t blocksize;	/* block size (in bytes) */
  size_t blockno;	/* number of blocks */
  size_t newblocksize;	/* new block size */
  size_t newblockno;	/* new block count */
  size_t rsize;		/* real size in bytes */
  void   *ptr;		/* the actual value */
} MEMSLOT;

/* struct MEMSLOT memory_slots[]
*    static array containing for each slot the actual blocksize, number
*    of blocks, and a pointer to the actual location. The location can
*    change when reallocating any other memory slot. */
static MEMSLOT memory_slots[M_MSLOTSTOTAL]; /* memory slots */

/* bool OUT_OF_MEMORY
*    flag indicating whether we are out of memory. */
#define OUT_OF_MEMORY	dd_stats.out_of_memory

/* void report_memory_usage(void)
*    for each used slot report the blocksize, number of blocks,
*    total memory used by this slot, and the actual pointer. */
void report_memory_usage(void)
{int i; MEMSLOT *ms;
    report(R_txt,"vertex and facet storage%s:\n",
        OUT_OF_MEMORY ? " (out of memory)":"");
    for(i=0,ms=&memory_slots[0];i<M_MSLOTSTOTAL;i++,ms++) if(ms->ptr) {
        report(R_txt,"M   slot %2d: blocksize=%6zu, no=%6zu, total=%9zu, ptr=%p\n",
           i,ms->blocksize,ms->blockno,ms->rsize,ms->ptr);
    }
}

/* void initialize_memory_slots(void)
*    clear all entries in memory_slots[] */
#define initialize_memory_slots()	\
    memset(&memory_slots[0],0,sizeof(memory_slots))

/* void init_main_slot(memslot_t slot, int nno, int nsize)
*    initializes a main slot by allocating the requested memory.
*     slot:    memory slot to be initialized
*     nno:     number of initial blocks
*     nsize:   initial block size
*/
static void init_main_slot(memslot_t slot, size_t nno, size_t nsize)
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
           "Allocmem: out of memory for slot=%d, blocksize=%zu, n=%zu\n",
           slot,nsize,nno);
        OUT_OF_MEMORY=1;
        return;
    }
    memset(ms->ptr,0,total);
    return;
}

/* void yalloc(type,slot,n,bsize)
*    initializes a main memory slot by requesting n blocks, where
*    each block is an array of bsize elements of the given type. */
#define yalloc(type,slot,n,bsize)	\
    init_main_slot(slot,n,(bsize)*sizeof(type))


/* void request_main_mem(memslot_t slot, int nno, int nsize)
*    records the requested block count and block size for a main
*    slot. The actual memory allocation is done by reallocmem() */
static inline void request_main_mem(memslot_t slot, size_t nno, size_t nsize)
{MEMSLOT *ms;
    ms=&memory_slots[slot];
    if(ms->blocksize==nsize && nno<=ms->blockno) return;
    ms->newblocksize=nsize;
    ms->newblockno=nno;
}    

/* void yrequest(type,slot,n,bsize)
*    request memory at a main slot; should be followed by calling
*    reallocmem(). */
#define yrequest(type,slot,n,bsize)	\
    request_main_mem(slot,n,(bsize)*sizeof(type))

/* int reallocmem(void)
*    allocate previously requested memory; if successful, adjust
*    blocks to the given count and size, and clear the new part.
*    Return 1 if out of memory, otherwise return 0. */
static int reallocmem(void)
{MEMSLOT *ms; int j,success; size_t total; void *ptr;
    for(j=0,success=1,ms=&memory_slots[0];
        success && j<M_MAINSLOTS; j++,ms++) if(ms->newblocksize){
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
        j<M_MAINSLOTS; j++,ms++) if(ms->newblocksize){
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

/* void* init_temp_slot(memslot_t slot, int nno, int nsize)
*    requests nno blocks, each of size nsize at the given slot.
*    The returned memory is not cleared; should check OUT_OF_MEMORY */
static void* init_temp_slot(memslot_t slot, size_t nno, size_t nsize)
{size_t total; MEMSLOT *ms;
    if(OUT_OF_MEMORY) return NULL;
    ms=&memory_slots[slot];
    ms->blocksize=nsize; ms->blockno=nno;
    total=nno*nsize;
    if(total <= ms->rsize) return ms->ptr;
    if(ms->ptr){ 
        free(ms->ptr);
        dd_stats.total_memory -= ms->rsize;
    }
    ms->rsize=total;
    dd_stats.total_memory += total;
    ms->ptr = malloc(total);
    if(!ms->ptr){
        report(R_fatal,"Out of memory for slot=%d, blocksize=%zu, n=%zu\n",
            slot,nsize,nno);
        OUT_OF_MEMORY=1;
        return NULL;
    }
    return ms->ptr;
}

#define talloc(type,slot,n,bsize)	\
    (type*) init_temp_slot(slot,n,(bsize)*sizeof(type))

/* void* request_temp_mem(memslot_t slot,size_t nno)
*    request more blocks for the initialized temporary memory slot */
static inline void* request_temp_mem(memslot_t slot,size_t nno)
{size_t total; MEMSLOT *ms; void *ptr;
    ms=&memory_slots[slot];
    if(OUT_OF_MEMORY) return ms->ptr;
    total = ms->blocksize*nno;
    if(total<=ms->rsize) return ms->ptr;
    dd_stats.memory_allocated_no++;
    ptr=realloc(ms->ptr,total);
    if(!ptr){ OUT_OF_MEMORY=1; return ms->ptr; }
    dd_stats.total_memory += total - ms->rsize;
    ms->rsize=total;
    ms->blockno=nno;
    ms->ptr=ptr;
    return ms->ptr;
}

#define trequest(type,slot,n)		\
    (type*) request_temp_mem(slot,n)

/* yalloc   => init main slot
   yrequest => change block count/size in main slot
   talloc   => init a temporary slot
   trequest => request memory for a temporary slot
*/

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
#define packsizelog	3	/* sizeof(BITMAP_t)== 1<<packsizelog */
#endif

#define packsize	(1<<packsizelog)	/* 4 or 8 bytes */
#define packshift	(packsizelog+3) 	/* in bits (5 or 6 )*/
#define packmask	((1<<packshift)-1)	/* 31 or 63 */
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
*    number of BITMAP_t words in a vertex bitmap */
static int
  VertexSize,		// size of the vertex coordinate block
  NextVertex,		// next free slot for a vertex
  MaxVertices,		// upper bound for vertex numbers
  VertexBitmapBlockSize;// size of the vertex bitmap block

/* double *StoreVertices
*    where vertices are stored; each vertex is a block of VertexSize
*    doubles, and there is a room for MaxVertices vertex */
#define StoreVertices	\
    ((double*)memory_slots[M_vertexcoord].ptr)

/* BITMAP_t *FacetAdjList
*    bitmap with VertexBitmapBlockSize blocksize; stores the
*    adjacency list of facets in MaxFacets blocks. */
#define FacetAdjList	\
    ((BITMAP_t*)memory_slots[M_facetadj].ptr)

/* BITMAP_t *FACET_adj(fno)
*    the vertex adjacency list of the facet fno < MaxFacets */
#define FACET_adj(fno)	\
    (FacetAdjList+((fno)*VertexBitmapBlockSize))

/* double *VERTEX_coord(vno)
*    The block of VertexSize coordinates of the vertex vno */
#define VERTEX_coord(vno)\
    (StoreVertices+((vno)*VertexSize))

/* BITMAP_t *VERTEX_adj(vno)
*    the block of the facet adjacency list of vertex vno */
#define VERTEX_adj(vno)	\
    (VertexAdjList+((vno)*FacetBitmapBlockSize))

/********************************************************************
* F A C E T S
*
* int FacetSize
*    number of double values to store a facet equation
* int NextFacet
*    living facets has index smaller than NextFacet
* int MaxFacets
*    number of bits in a facet bitmap, NextFacet <= MaxFacets
* int FacetBitmapBlockSize
*    number of BITMAP_t words in a facet bitmap */
static int
  FacetSize,		// size of the facet equation block
  NextFacet,		// next free slot for a facet
  MaxFacets,		// number of bits in a facet bitmap
  FacetBitmapBlockSize;	// size of the facet bitmap block

/* double *StoreFacets
*    where facets are stored as block of FacetSize doubles. There is
*    space to store MaxFacets such blocks. */
#define StoreFacets	\
    ((double*)memory_slots[M_facetcoord].ptr)

/* BITMAP_t *VertexAdjList, *FacetLiving, *FacetFinal
*    bitmaps with FacetBitmapBlockSize BITMAP_t blocksize. VertexAdjList
*    stores the adjacency list of vertices; there are MaxVertices blocks.
*    Others are facet attributes and they occupy a single block. Each block
*    can store up to MaxFacets bits. FacetLiving tells whether a facet is 
*    among the facets of the current approximation; FacetFinal tells
*     whether the facet belongs to the final polytope. */
#define VertexAdjList	\
    ((BITMAP_t*)memory_slots[M_vertexadj].ptr)
#define FacetLiving	\
    ((BITMAP_t*)memory_slots[M_facetliving].ptr)
#define FacetFinal	\
    ((BITMAP_t*)memory_slots[M_facetfinal].ptr)

/* double *FACET_coords(fno)
*    the facet equation; each equation occupies FacetSize doubles,
*    and fno < MaxFacets */
#define FACET_coords(fno)\
    (StoreFacets+((fno)*FacetSize))

/***********************************************************************
* Some auxiliary bitmap functions
*
* bool is_livingFacet(fno)
*     whether facet fno is an actual facet */
#define is_livingFacet(idx)	\
	((FacetLiving[(idx)>>packshift]>>((idx)&packmask))&1)

/* bool is_finalFacet(fno)
*     whether facet fno is a facet of the final polytope */
#define is_finalFacet(idx)	\
	((FacetFinal[(idx)>>packshift]>>((idx)&packmask))&1)

/* void clear_facet_in_Fliving(fno)
*    clear bit fno in the bitmap FacetLiving */
#define clear_facet_in_Fliving(fno) \
    FacetLiving[(fno)>>packshift] &= ~(BITMAP1<<((fno)&packmask))

/* void set_facet_in_Fliving(fno)
*    set bit fno in bitmap FacetLiving */
#define set_facet_in_Fliving(fno) \
    FacetLiving[(fno)>>packshift] |= BITMAP1<<((fno)&packmask)

/* void set_facet_in_Ffinal(fno)
*    set bit fno in bitmap Facetfinal */
#define set_facet_in_Ffinal(fno)  \
    FacetFinal[(fno)>>packshift] |= BITMAP1<<((fno)&packmask)

/* clear_FacetAdj_all(facetno)
*    clear the adjacency list of this facet */
#define clear_FacetAdj_all(fno)	  \
    memset(FACET_adj(fno),0,VertexBitmapBlockSize*sizeof(BITMAP_t))

/* set_FacetAdj(facetno,vertexno)
*    set bit vertexno in the adjacency list of facet facetno */
#define set_FacetAdj(fno,vno)	  \
    FACET_adj(fno)[(vno)>>packshift] |= BITMAP1 << ((vno)&packmask)

/* clear_VertexAdj_all(vertexno)
*    clear the adjacency list of this vertex */
#define clear_VertexAdj_all(vno)  \
    memset(VERTEX_adj(vno),0,FacetBitmapBlockSize*sizeof(BITMAP_t))

/* set_VertexAdj(vertexno,facetno)
*    set bit facetno in the adjacency list of vertex vertexno */
#define set_VertexAdj(vno,fno)	  \
    VERTEX_adj(vno)[(fno)>>packshift] |= BITMAP1 << ((fno)&packmask)

/* intersect_VertexAdj_Fliving(vertexno)
*    clear bits in the adjacency list of vertexno which are not in
*    FacetLiving */
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
*
* int bitcnt_mask, bitcnt_shift
*    number of bits stored for the last bitcnt_shift many bits
*
* char bitcnt[0..bitcnt_mask]
*    number of bits in the index. Used to speed up counting bits
*    in a bitmap.
*
* int add_bitcount(BITMAP_t v, int *total)
*    add the number of bits set in v to *total
*/
DD_STATS dd_stats;

#define bitcnt_shift	10
#define bitcnt_mask	((1u<<bitcnt_shift)-1u)

/* how many bits are in 0..bitcnt_mask */
static char bitcnt[] = {
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10
};

#define add_bitcount(v,total)	\
    while(v){ total += bitcnt[(v)&bitcnt_mask]; (v)>>=bitcnt_shift; }

void get_dd_facetno(void) // get the number living and final facets
{int i; BITMAP_t v;
    dd_stats.living_facets_no=-1;
    dd_stats.final_facets_no=-1;
    for(i=0;i<FacetBitmapBlockSize;i++){
        v=FacetLiving[i];
        if(v==~BITMAP0){ dd_stats.living_facets_no+=(1<<packshift); }
        else add_bitcount(v,dd_stats.living_facets_no);
        v=FacetFinal[i];
        if(v==~BITMAP0){ dd_stats.final_facets_no +=(1<<packshift); }
        else add_bitcount(v,dd_stats.final_facets_no);
    }
    if(dd_stats.final_facets_no<0) dd_stats.final_facets_no=0;
    if(dd_stats.living_facets_no<0) dd_stats.living_facets_no=0;
}
/***********************************************************************
* Interrogate facets and vertices
*
* int mrandom(v)
*    return a random integer between 0 and v-1 approximately uniformly */
static inline int mrandom(int v)
{ return v<=1?0 : (random()&0x3fffffff)%v; }

/* int get_next_facet(from)
*    return the next living but not final facet number starting at
*    from. If no such vertices are, return -1.
*    When from=-1 and RandomFacet is set, pick the starting number
*    randomly. */
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

/* void mark_facet_as_final(fno)
*    mark the facets as final, i.e. facet of the final polytope. */
void mark_facet_as_final(int fno)
{    FacetFinal[fno>>packshift] |= BITMAP1 << (fno&packmask); }

/* int vertex_num(void)
*    number of vertices added so far (except the ideal ones) */
int vertex_num(void)
{ return NextVertex-DIM; }

/* int facet_num(void)
*    return the number of the facets of the actual approximating polytope;
*    this is the number of facets of the final polytope when the algorithm
*    terminates */
int facet_num(void)
{int i,j;
    j=0; for(i=1;i<NextFacet;i++)if(is_livingFacet(i)) j++;
    return j;
}

/* double round(double x)
*    round x to an integer of closer than SCALE_EPS */
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

/* bool closetoint(double x)
*    check if the argument is close to a (small) integer value */
inline static int closetoint(double x)
{   if(x<0.0) x=-x;
    x -= (int)(x+0.5); 
    return (x<1e-8 && x>-1e-8);
}

/* int gcd(int a, int b)
*  int lcm(int a, int b)
*    compute the gcd and lcm of two integers */
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

/* int denum(double x)
*    check if x is close to a rational with denominator at most 1000.
*    return the value of the denominator */
static int denum(double x)
{int i;
    for(i=1;i<1000;i++)if(closetoint(i*x)) return i;
    return 1;
}

/* void get_facet_into(fno, double v[0:dim])
*    store the facet equation to the provided space. The equation is 
*    scaled and rounded to the closet integer when the difference is 
*    small */
void get_facet_into(int fno, double *v)
{int j,d;
    d=1; for(j=0;d<130000 && j<=DIM;j++)d=lcm(d,denum(FACET_coords(fno)[j]));
    for(j=0;j<=DIM;j++)v[j]=round(d*FACET_coords(fno)[j]);
}

/* char *formatvalue(double v)
*    formatted output of a vertex coordinate */
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

/* void print_vertex(report_type, double dir, double coords[0:DIM-1](
*    report the vertex using fractional or floating format. */
void print_vertex(report_type where, double dir, const double v[/*0:DIM-1*/])
{int j;
    for(j=0;j<PARAMS(ProblemObjects);j++){ // DIM might not be set
        report(where," %s",formatvalue(dir*v[j]));
    }
    report(where,"\n");
}

/* void print_vertices(report_type where)
*    report all stored vertices using print_vertex(). */
void print_vertices(report_type where)
{int i; double dir;
    dir = PARAMS(Direction) ? -1.0 : +1.0;
    for(i=DIM;i<NextVertex;i++){
        report(where,"V ");
        print_vertex(where,dir,VERTEX_coord(i));
    }
}

/* void print_facets(report_type where)
*    report all facets using scaling and rounding format. */
void print_facets(report_type where)
{int i,j,d;
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
    MaxFacets = DD_INITIAL_FACETNO;
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
    // allocate memory for the main slots
    initialize_memory_slots();
    yalloc(double,M_vertexcoord,MaxVertices,VertexSize); // StoreVertices
    yalloc(double,M_facetcoord,MaxFacets,FacetSize); // StoreFacets
    yalloc(BITMAP_t,M_vertexadj,MaxVertices,FacetBitmapBlockSize); // VertexAdjList
    yalloc(BITMAP_t,M_facetadj,MaxFacets,VertexBitmapBlockSize); // FacetAdjList
    yalloc(BITMAP_t,M_facetliving,1,FacetBitmapBlockSize); // FacetLiving
    yalloc(BITMAP_t,M_facetfinal,1,FacetBitmapBlockSize); // FacetFinal

//    yalloc(double,M_facetdist,MaxFacets,1); // FACET_dist
//    yalloc(int,M_facetposneg,MaxFacets,1); // FacetPostnegList
//    yalloc(int,M_vertexlist,MaxVertices,1); // VertexList
//    yalloc(BITMAP_t,M_facetwork,1,FacetBitmapBlockSize); // FacetWork
//    yalloc(double,M_vertexarray,VERTEXARRAY_STEPSIZE,DIM+1); // VertexArray

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
*    add space for DD_VERTEX_ADDBLOCK<<packshift more vertices, then
*    expand StoreVertices, VertexAdjList;
*    adjust VertexBitmapBlockSize by DD_VERTEX_ADDBLOCK
*    Called at the very beginning of a new iteration when necessary. */
static void allocate_vertex_block(void)
{   // extend Vertices and FacetBitmap to accommodate more stuff
    dd_stats.vertices_allocated_no ++;
    dd_stats.vertices_allocated += (DD_VERTEX_ADDBLOCK<<packshift);
    MaxVertices += (DD_VERTEX_ADDBLOCK<<packshift);
    VertexBitmapBlockSize += DD_VERTEX_ADDBLOCK;
    // tell the memory handling part how much storage space we would need.
    yrequest(double,M_vertexcoord,MaxVertices,VertexSize);
    yrequest(BITMAP_t,M_vertexadj,MaxVertices,FacetBitmapBlockSize);
//    yrequest(int,M_vertexlist,MaxVertices,1);
    yrequest(BITMAP_t,M_facetadj,MaxFacets,VertexBitmapBlockSize);
    // and do the reallocation
    if(reallocmem()){  // out of memory, don't increase the values
        MaxVertices -= (DD_VERTEX_ADDBLOCK<<packshift);
        VertexBitmapBlockSize -= DD_VERTEX_ADDBLOCK;
    }
}

/* void allocate_facet_block(int count)
*    add space for count more facets. Increase storage place in 
*    StoreFacets and FacetAdjList, add more bits to facet bitmaps */
static void allocate_facet_block(int count)
{int total;
    if(OUT_OF_MEMORY) return;
    for(total=DD_FACET_ADDBLOCK<<packshift;total<count; total += DD_FACET_ADDBLOCK<<packshift);
    dd_stats.facets_allocated_no ++;
    dd_stats.facets_allocated += total;
    MaxFacets += total;
    FacetBitmapBlockSize = (MaxFacets+packmask)>>packshift;
    yrequest(double,M_facetcoord,MaxFacets,FacetSize);
    yrequest(BITMAP_t,M_vertexadj,MaxVertices,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_facetadj,MaxFacets,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_facetliving,1,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_facetfinal,1,FacetBitmapBlockSize);
    if(reallocmem()){ // out of memory
        FacetBitmapBlockSize -= (MaxFacets+packmask)>>packshift;
        MaxFacets -= total;
    }
}

/***********************************************************************
* Compute facet equation from the vertices it is adjacent to.
*
* bool solve_lineq(int d, double VertexArray[dim+1,d])
*    solve the system of homogeneous linear equation stored in
*    VertexArray[dim+1,d]. Use Gauss elimination: get the largest value
*    in the first column A[0,], swap this row and the first row; and
*    subtract from all rows so that the first column will be all zero,
*    etc. The matrix rank should be exactly dim (so that we have a
*    non-trivial solution). The solution is returned in A[0..dim].
*    PARAMS(LineqEps) is the rank threshold.
*    Returns 1 if the matrix rank is not dim; otherwise returns zero. */
static int solve_lineq(int d, double *VertexArray)
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

/* void recalculate_facet_eq(int fno, BITMAP_t *adj, double *coords, memslot_t tmp)
*    calculate facet equation form the vertices adjacent to it. Complain
*    if out of memory, the system is degenerate, or the old and new coeffs
*    are too far from each other. */
static void recalculate_facet_eq(int fno,BITMAP_t *facetadj, double *facetcoords,
    memslot_t MemTmp)
{double v,s; int i,j,vno; int amax,an; BITMAP_t fc;
 double *VertexArray;    
#define A(i,j)	VertexArray[(i)*(DIM+1)+(j)]
    // collect all vertices adjacent to facet fno
    amax=VERTEXARRAY_STEPSIZE; /* we have at least that many slots */
    VertexArray = talloc(double,MemTmp,amax,DIM+1);
    an=0; vno=0;
    for(i=0;i<VertexBitmapBlockSize;i++){
        j=vno; fc=facetadj[i];
            while(fc){
            if(fc&1){ /* store vertex j to A */
                if(an>=amax){
                    amax+=VERTEXARRAY_STEPSIZE;
                    VertexArray = trequest(double,MemTmp,amax);
                    if(OUT_OF_MEMORY) return; // out of memory
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
    if(solve_lineq(an,VertexArray)){
        report(R_err,"recalculate_facet: adjacency vertex list of facet %d is degenerate\n",fno);
        dd_stats.numerical_error++;
        return;
    }
    s=0.0; for(i=0;i<DIM;i++) s+= VertexArray[i]; s=1.0/s;
    for(i=0;i<=DIM;i++){
        VertexArray[i]*=s;
        v=facetcoords[i]-VertexArray[i];
        if(v>PARAMS(FacetRecalcEps) || v < -PARAMS(FacetRecalcEps)){
            report(R_warn,"recalculate_facet: numerical instability at facet %d, coord %d (%lg)\n", fno,i,v);
            dd_stats.instability_warning++;
        }
        facetcoords[i]=VertexArray[i];
    }
#undef A
}

/* void recalculate_facets(void)
*    go over all facets, and recalculate their equations. */
void recalculate_facets(void)
{int fno; 
    for(fno=1;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        recalculate_facet_eq(fno,FACET_adj(fno),FACET_coords(fno),M_vertexarray);
    }
}

/***********************************************************************
* double vertex_distance(double *coords, int fno)
*    compute the distance of the vertex given as the first argument
*    from the facet given with its number */
#define DD_EPS_EQ	PARAMS(PolytopeEps)

inline static double vertex_distance(double *coords, int fno)
{double d=0.0; int i; double *fcoords;
    fcoords=FACET_coords(fno);
    for(i=0;i<DIM;i++){
       d+= (*coords)*(*fcoords);
       coords++; fcoords++;
    }
    d += (*fcoords);
    return d;
}

/* bool store_vertex(double *coords)
*    check if the vertex is new, and if yes, store it, and disable
*    further vertex processing. Return 1 if this vertex was not seen
*    before; and zero otherwise. */
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

/* int probe_vertex(double *coords)
*    return the number of negative facets. */
int probe_vertex(double *coords)
{int fno,negfacets;
    dd_stats.probevertex++;
    negfacets=0;
    for(fno=0;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        if(vertex_distance(coords,fno) < -DD_EPS_EQ) {
            negfacets++;
        }
    }
    return negfacets;
}

/***********************************************************************
* Routines for vertex enumeration
*
* An iteration means adding a new vertex to the existing polytope.
* The main routine is add_new_vertex(). It allocates memory to 
* StoreFacetDists and FacetPosnegList, used in the outermost loop.
* Then allocates memory for FacetWork and VertexList, both used
* in the innermost loop but does not change their size. Finally
* set the memory blocks for the newly created facets.
*
* V A R I A B L E S 
*
* double StoreFacetDists[MaxFacets]
*    for each facet it contains the distance of the facet and
*    the new vertex. It is fixed during the iteration.
*
* int *FacetPosnegList[MaxFacets]
*    indices of positive (starting from the beginning) and negative
*    (starting from the end) facets when compared to the newly
*    added vertex. It is fixed during the iteration. */

static double* StoreFacetDists;
static int*    FacetPosnegList;
/* double FACET_dist(fno)
*    the distance associated with facet fno < MaxFacets */
#define FACET_dist(fno)	\
    StoreFacetDists[fno]

/* BITMAP_t* FacetWork
*    auxiliary facet bitmap, it stores living facets minus the two
*    facets whose intersection we are testing. Used in is_ridge()
*    only. The size is fixed, the content changes.
*
* int VertexList[MaxVertices]
*    list to store vertex indices in the intersection of the two
*    facets. Used in is_ridge() only. The size is fixed, the
*    content changes. */
static BITMAP_t* FacetWork;
static int* VertexList;

/*
* int NewFacet
*    number of newly created facets
* int MaxNewFacets
*    have storage space for that many newly created facets */
static int 
  NewFacet,			// number of newly created facets
  MaxNewFacets;			// available space

/* double NewStoreFacets
*    coordinates of the newly created facets are stored here.
*    There is a space to store MaxNewFacets such blocks. */
static double* NewStoreFacets;

/* BITMAP_t NewFacetAdjList
*    the vertex adjacency list of the newly generated facet */
static BITMAP_t* NewFacetAdjList;

/* double *NewFACET_coords(fno)
*    the facet equation; each equation occupies FacetSize doubles,
*    and fno < MaxFacets */
#define NewFACET_coords(fno)\
    (NewStoreFacets+((fno)*FacetSize))

/* BITMAP_t *NewFACET_adj(fno)
*    the vertex adjacency list of the facet fno < MaxFacets */
#define NewFACET_adj(fno)	\
    (NewFacetAdjList+((fno)*VertexBitmapBlockSize))

/* set_NewFACET_adj(facetno,vertexno)
*    set bit vertexno in the adjacency list of the new facet facetno */
#define set_NewFACET_adj(fno,vno)	  \
    NewFACET_adj(fno)[(vno)>>packshift] |= BITMAP1 << ((vno)&packmask)

/***********************************************************************
* Some bitmap manipulating routines
*
* void copy_Fliving_to_Fwork()
*    copy the bitmap FacetLiving to FacetWork */
#define copy_Fliving_to_Fwork()	  \
    memcpy(FacetWork,FacetLiving,FacetBitmapBlockSize*sizeof(BITMAP_t))

/* void clear_facet_in_Fwork(fno)
*    clear bit fno in the bitmap FacetWork */
#define clear_facet_in_Fwork(fno) \
    FacetWork[(fno)>>packshift] &= ~(BITMAP1<<((fno)&packmask))

/***********************************************************************
* void move_newfacet_to(double *coords, BITMAP_t *adj, fno)
*   copy the facet at NextFacet to the given place, and decrease
*   NextFacet */
#define move_newfacet_to(coords,adj,fno)	\
  { memcpy(FACET_coords(fno),coords,FacetSize*sizeof(double)); \
    memcpy(FACET_adj(fno),adj,VertexBitmapBlockSize*sizeof(BITMAP_t)); }

/***********************************************************************
* Get the next facet number
* int get_new_facetno()
*    return the new facet number, asking for memory if necessary.
*    The return value is -1 if cannot allocate more memory.
*/
inline static int get_new_facetno(void)
{int i;
    i=NewFacet; 
    if(NewFacet>=MaxNewFacets){ // no more space, ask memory
        MaxNewFacets += DD_FACET_ADDBLOCK<<packshift;
        NewStoreFacets = trequest(double,M_newfacetcoord,MaxNewFacets);
        NewFacetAdjList = trequest(BITMAP_t,M_newfacetadj,MaxNewFacets);
        if(OUT_OF_MEMORY) {
            MaxNewFacets -= DD_FACET_ADDBLOCK<<packshift;
            return -1;
        }
    }
    NewFacet++; return i;
}

/***********************************************************************
* int facet_intersection(f1,f2)
*    intersect the vertex adjacency lists of f1 and f2; return the
*    number of vertices adjacent to both f1 and f2 */
inline static int facet_intersection(int f1, int f2)
{int i,total; register BITMAP_t *L1,*L2, v;
    total=0; L1=FACET_adj(f1); L2=FACET_adj(f2);
    for(i=0;i<VertexBitmapBlockSize;i++,L1++,L2++){
        v=(*L1)&(*L2);
        add_bitcount(v,total);
    }
    return total;
}

/* int is_ridge(f1,f2)
*    combinatorial test to check whether the intersection of f1 and
*    f2 forms a ridge. If there are <DIM-1 vertices both in f1 and
*    f2 then it is not a ridge. Otherwise store vertex indices in
*    VertexList, then take the intersection of living facets and 
*    the facet adjacency lists in VertexList. This is a ridge if 
*    and only if the intersection contains f1 and f2 only. */
inline static int is_ridge(int f1,int f2)
{int vertexno,i,j,vlistlen; BITMAP_t v;
    if(facet_intersection(f1,f2) < DIM-1){
         return 0; // no - happens oftern
    }
    copy_Fliving_to_Fwork();
    // only f1 and f2 should remain, if any other remains, the answer is no
    clear_facet_in_Fwork(f1); clear_facet_in_Fwork(f2);
    // loop through the intersection of FACET_adj(f1) and FACET_adj(f2)
    // store these vertices in VertexList
    vertexno=0; vlistlen=0;
    for(i=0;i<VertexBitmapBlockSize;i++){
        if((v=FACET_adj(f1)[i] & FACET_adj(f2)[i])){
            j=vertexno;
            while(v){
                while((v&7)==0){ j+=3; v>>=3; }
                if(v&1){ VertexList[vlistlen]=j; vlistlen++; }
                j++; v>>=1;
            }
        }
        vertexno += (1<<packshift);
    }
    /* this code works only when vlistlen>=2. This holds for sure
       when DIM >=3 */
    for(i=0;i<FacetBitmapBlockSize;i++) if(
      (v=FacetWork[i] & VERTEX_adj(VertexList[0])[i])
      && (v &= VERTEX_adj(VertexList[1])[i])){
        for(j=2;j<vlistlen;j++){
            v &= VERTEX_adj(VertexList[j])[i];
        }
        if(v) return 0; // no
    }
    return 1; // yes
}

/* void create_new_facet(f1,f2,vno)
*    create a new facet through the ridge formed by f1 and f2, and
*    the given vertex vno. f1*v is negative, f2*v is positive, and these
*    values are stored in FACET_dist[]. Recalculate the equation when
*    ExactFacetEq parameter is set */
inline static void create_new_facet(int f1, int f2, int vno)
{int newf; double d1,d2,d; int i;
    newf=get_new_facetno();
    if(newf<0) return; // no memory
    for(i=0;i<VertexBitmapBlockSize;i++)
        NewFACET_adj(newf)[i] = FACET_adj(f1)[i] & FACET_adj(f2)[i];
    set_NewFACET_adj(newf,vno);
    // compute the coefficients, f1<0, f2 >0
    if(f2==0){ // ideal facet
        for(i=0;i<=DIM;i++){
            NewFACET_coords(newf)[i]=FACET_coords(f1)[i];
        }
        NewFACET_coords(newf)[DIM] -= FACET_dist(f1);
    } else {   // f1(v)=-d1, f2(v)=d2; (d2*f1 + d1*f2)(v)=0.0
        d1=-FACET_dist(f1); d2=FACET_dist(f2);
        d=1.0/(d1+d2); d1 *= d; d2 *= d;
        for(i=0;i<=DIM;i++){
            NewFACET_coords(newf)[i]=d2*FACET_coords(f1)[i]+d1*FACET_coords(f2)[i];
        }
    }
    if(PARAMS(ExactFacetEq))
       recalculate_facet_eq(MaxFacets+newf,NewFACET_adj(newf),NewFACET_coords(newf),M_vertexarray);
}

/* void search_ridges_with(f1,vno)
*    Vertex vno is on the negative side of f1; go over all positive
*    facets and call create_new_facet() for ridges */
inline static void search_ridges_with(int f1, int newvertex)
{int *f2,j;
    for(f2=FacetPosnegList; (j=*f2)>=0; f2++)
        if(is_ridge(f1,j)) create_new_facet(f1,j,newvertex);
}

/* void search_ridges_DIM2(f1,vno)
*    same as search_ridges_with() but for the case DIM==2 */
inline static void search_ridges_DIM2(int f1, int newvertex)
{int *f2,j; int i,total;  BITMAP_t *L1, *L2, v;
    for(f2=FacetPosnegList; (j=*f2)>=0; f2++){
        total=0; L1=FACET_adj(f1); L2=FACET_adj(j);
        for(i=0;i<VertexBitmapBlockSize;i++,L1++,L2++){
            v=(*L1)&(*L2);
            add_bitcount(v,total);
        }
        // facets f1 and j intersect in a vertex
        if(total!=0) create_new_facet(f1,j,newvertex);
    }
}

/* void make_facet_living(fno)
*    set the "living" flag for this facet, and add it to the adjacency
*    list of all vertices it is adjacent to. */
static void make_facet_living(int fno)
{int vno,i,j; BITMAP_t fc;
    set_facet_in_Fliving(fno);
    vno=0; for(i=0;i<VertexBitmapBlockSize;i++){
        j=vno; fc=FACET_adj(fno)[i];
        while(fc){
            while((fc&7)==0){ j+=3; fc>>=3; }
            if(fc&1){set_VertexAdj(j,fno);}
            j++; fc>>=1;
        }
        vno+=(1<<packshift);
    }
}

/* void add_new_vertex(double vertex[0:dim-1])
*    add a new vertex which is outside the convex hull of the present
*    approximation. Split facets into positive, negative and zero sets
*    depending where the new vertex is. For each pair of positive/
*    negative facets, check if it is a ridge; if yes, add a new facet.
*    Then throw away negative facets, and add the newly created facets
*    to the approximation. */

/** add a new vertex to the approximation **/
void add_new_vertex(double *coords)
{int i,j,fno; int thisvertex; BITMAP_t fc; double d;
 int *PosIdx, *NegIdx;
    dd_stats.iterations++;
    if(NextVertex >= MaxVertices){
        allocate_vertex_block();
    }
    // memory blocks for the outer loop
    StoreFacetDists = talloc(double,M_facetdist,MaxFacets,1);
    FacetPosnegList = talloc(int,M_facetposneg,MaxFacets,1);
    // blocks for the inner loop
    FacetWork = talloc(BITMAP_t,M_facetwork,1,FacetBitmapBlockSize);
    VertexList = talloc(int,M_vertexlist,MaxVertices,1);
    // space for the new facets
    NewFacet=0; MaxNewFacets=DD_INITIAL_FACETNO;
    NewStoreFacets = talloc(double,M_newfacetcoord,MaxNewFacets,FacetSize);
    NewFacetAdjList = talloc(BITMAP_t,M_newfacetadj,MaxNewFacets,VertexBitmapBlockSize);
    if(OUT_OF_MEMORY) return;

    thisvertex=NextVertex; NextVertex++;
    for(i=0;i<DIM;i++){ VERTEX_coord(thisvertex)[i]=coords[i]; }
    clear_VertexAdj_all(thisvertex); // clear the adjacency list
    // split facets into positive, negative and zero parts
    dd_stats.facet_pos=0; dd_stats.facet_neg=0; dd_stats.facet_zero=0;
    PosIdx = FacetPosnegList; // this goes ahead
    NegIdx = FacetPosnegList+MaxFacets; // this goes backward
    for(fno=0;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        d=FACET_dist(fno)=vertex_distance(coords,fno);
        if(d>DD_EPS_EQ){ // positive side 
            *PosIdx = fno; ++PosIdx;
            dd_stats.facet_pos++;
        } else if(d<-DD_EPS_EQ){ // negative side
            --NegIdx; *NegIdx=fno;
            dd_stats.facet_neg++;
        } else { // this is adjacent to our new vertex
            set_VertexAdj(thisvertex,fno);
            set_FacetAdj(fno,thisvertex);
            dd_stats.facet_zero++;
        }
    }
    *PosIdx = -1; --NegIdx; *NegIdx=-1;
    // add statistics:
    d=(double)dd_stats.facet_pos*(double)dd_stats.facet_neg;
    if(dd_stats.max_tests<d)dd_stats.max_tests=d;
    dd_stats.avg_tests = (dd_stats.iterations-1)*dd_stats.avg_tests*(1.0/dd_stats.iterations)
      + (1.0/dd_stats.iterations)*d;
    // this is the critical part: for each positive/negative facet pair
    // compute whether this is a new facet.
    // go through all negative facets
    NegIdx = FacetPosnegList+MaxFacets;
    if(DIM<=2){ // search_ridges_with() assumes DIM>=3 
        while((j=*--NegIdx)>=0) search_ridges_DIM2(j,thisvertex);
    } else {
        while((j=*--NegIdx)>=0) search_ridges_with(j,thisvertex);
    }
    if(OUT_OF_MEMORY || dobreak){
         dd_stats.facet_new=0;
         return;
    }
    // new facets are 0 .. NewFacet in NewFACET()
    // calculate the number of facets of the new polytope
    dd_stats.facet_new=NewFacet;
    i = dd_stats.facet_zero + dd_stats.facet_pos + dd_stats.facet_new;
    if(dd_stats.max_facets < i) dd_stats.max_facets = i;
    if(dd_stats.max_facetsadded<dd_stats.facet_new){
        dd_stats.max_facetsadded=dd_stats.facet_new; 
    }
    dd_stats.avg_facetsadded = (dd_stats.iterations-1)*dd_stats.avg_facetsadded*(1.0/dd_stats.iterations)
         + (1.0/dd_stats.iterations)*dd_stats.facet_new;
    // delete negative facets from FacetLiving
    NegIdx = FacetPosnegList+MaxFacets;
    while((j=*--NegIdx)>=0) clear_facet_in_Fliving(j);
    // move new facets to NextFacet; this part can be skipped 
    while(NewFacet>0 && NextFacet<MaxFacets){
        NewFacet--;
        move_newfacet_to(NewFACET_coords(NewFacet),NewFACET_adj(NewFacet),NextFacet);
        make_facet_living(NextFacet);
        NextFacet++;
    }
    if(NewFacet==0) return;
    // NextFacet == MaxFacets, no more space here
    dd_stats.facets_compressed_no ++;
    // clear adjacency list of vertices
    for(i=0;i<NextVertex;i++){ // clear the adjacency list of vertices
        intersect_VertexAdj_Fliving(i);
    }
    // move new facets into free facet slots
    fno=0; for(i=0;i<FacetBitmapBlockSize;i++){
        j=fno; fc=~FacetLiving[i]; // complement ...
        while(fc){
            while((fc&7)==0){ j+=3; fc>>=3; }
            if((fc&1)){
                if(j>=MaxFacets) goto finish_compress;
                NewFacet--;
                move_newfacet_to(NewFACET_coords(NewFacet),NewFACET_adj(NewFacet),j); 
                make_facet_living(j);
                if(NewFacet==0) goto finish_compress;
            }
            j++; fc>>=1;
        }
        fno += (1<<packshift);
    }
  finish_compress:
    if(NewFacet>0){ // all is full, ask more space
        allocate_facet_block(NewFacet);
        // if no memory, throw the newly generated facets
        if(OUT_OF_MEMORY) NewFacet=0;
        while(NewFacet>0){
            // ASSERT(NextFacet < MaxFacet);
            NewFacet--;
            move_newfacet_to(NewFACET_coords(NewFacet),NewFACET_adj(NewFacet),NextFacet);
            make_facet_living(NextFacet);
            NextFacet++;
        }
        return;
    }
  fill_holes:
    while( fno<NextFacet && is_livingFacet(fno) ) fno++;
    while( NextFacet > fno && !is_livingFacet(NextFacet-1) ) NextFacet--;
    if(fno < NextFacet){
        NextFacet--;
        move_newfacet_to(FACET_coords(NextFacet),FACET_adj(NextFacet),fno);
        make_facet_living(fno);
        clear_facet_in_Fliving(NextFacet);
        if(is_finalFacet(NextFacet)) set_facet_in_Ffinal(fno);
        fno++; goto fill_holes;
    }
    // clear FacetFinal
    for(i=0;i<FacetBitmapBlockSize;i++) FacetFinal[i] &= FacetLiving[i];
    // clear the adjacency list of vertices
    for(i=0;i<NextVertex;i++){ // clear the adjacency list of vertices
        intersect_VertexAdj_Fliving(i);
    }
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

