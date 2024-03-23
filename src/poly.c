/** poly.c  --  combinatorial part using double description method **/

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

#ifdef USETHREADS
/************************************************************************
* Threads
*
* Creating new facets is split among PARAMS(Threads) threads. The 
* thread logic is the following:
*   MAIN THREAD:
*        loop {
*            ...
*            distribute work (or ask threads to terminate)
*            -- ThreadBarrierForking--
*            do 1/n of the parallel work
*            -- ThreadBarrierJoining --
*            parallel work is done
*            ...
*        }
*  EXTRA THREADS:
*        loop {
*            ...
*            -- ThreadBarrierForking --
*            optionally exit if instructed
*            get the work and do it
*            -- ThreadBarrierJoining --
*        }
*/

#include <pthread.h>
#include <sys/sysinfo.h>  /* get_nprocs() */

/************************************************************************
*
*       T H R E A D S
*
*************************************************************************
*
* Variables used by threads
*
* thread_data_t ThreadData[MAX_THREADS]
*    data block containing data on a single thread */
typedef struct {
    int		id;	/* thread number */
    int		quit;	/* set to 1 if the thread should stop */
    pthread_t	obj;	/* the thread itself */
} thread_data_t;

thread_data_t ThreadData[MAX_THREADS]; // 0 is unused

/* thread job, the argument is zero for the calling thread */
typedef void ThreadJob_t(int threadId);
static ThreadJob_t *ThreadJob;
#define ThreadNo        PARAMS(Threads)

/* pthread_barrier_t ThreadBarrierForking, ThreadBarrierJoining
*    barriers for synchronizing work */
static pthread_barrier_t ThreadBarrierForking;
static pthread_barrier_t ThreadBarrierJoining;

/* void *extra_thread(void *arg)
*    the code for extra threads - forward declaration */
static void *extra_thread(void *arg);

/* int create_threads(void)
*    create barriers, mutex, threads, and start them. Threads will
*    wait until they get their first assignment. */
int create_threads(void)
{int i,rc;
    if(PARAMS(ProblemObjects)<2){
        PARAMS(Threads)=1; return 0; // do not use threads in this case
    }
    if(PARAMS(Threads)<=0){ // figure out the number of CPU's
        PARAMS(Threads)=get_nprocs();
    }
    if(PARAMS(Threads)<=1){
        PARAMS(Threads)=1;
        return 0; // nothing to do
    }
    if(PARAMS(Threads)>MAX_THREADS){ PARAMS(Threads)=MAX_THREADS; }
    // barriers
    if((rc = pthread_barrier_init(&ThreadBarrierJoining, NULL, PARAMS(Threads)))){
        report(R_fatal,"error: pthread_barrier_init, rc: %d\n", rc);
        return 1;
    }
    if((rc = pthread_barrier_init(&ThreadBarrierForking, NULL, PARAMS(Threads)))){
        report(R_fatal,"error: pthread_barrier_init, rc: %d\n", rc);
        return 1;
    }
    // create and start threads, they will stop arriving at the barriers
    report(R_info,"Using %d threads\n", PARAMS(Threads));
    ThreadData[0].id=0;
    for(i=1;i<PARAMS(Threads);i++){
        ThreadData[i].id = i;
        ThreadData[i].quit = 0;
        if((rc = pthread_create(&(ThreadData[i].obj),NULL,extra_thread,&ThreadData[i]))){
            report(R_fatal,"error creating thread %d, rc: %d\n",i,rc);
            return 1;
        }
    }
    return 0;
}
/* void thread_execute(ThreadJob_t job)
*    execute "job" using all avaiable threads */
static void thread_execute(ThreadJob_t job)
{   if(PARAMS(Threads)<2){ // single thread
        job(0); return;
    }
    ThreadJob=job;
    pthread_barrier_wait(&ThreadBarrierForking); // start threads
      job(0); // main thread
    pthread_barrier_wait(&ThreadBarrierJoining); // wait until others finish
}

/* void *extra_thread(void *arg)
*    each thread waits, then calls the routine in ThreadJob, then sleeps */
static void *extra_thread(void *arg)
{thread_data_t *data = (thread_data_t*)arg;
int myId=data->id;
    report(R_info,"Thread %d started ...\n",myId);
    while(1){
       pthread_barrier_wait(&ThreadBarrierForking);
       if(data->quit) break; // stop
       ThreadJob(myId);
       pthread_barrier_wait(&ThreadBarrierJoining);
    }
    report(R_info,"Thread %d stopped\n",myId);
    return NULL;
}
/* void stop_threads(void)
*    tell extra threads to exit and wait for (reap) them.
*    This blocks until all threads exit. */
void stop_threads(void)
{int i;
    if(PARAMS(Threads)<2) return;
    report(R_info,"Stopping threads ...\n");
    for(i=1;i<PARAMS(Threads);i++) ThreadData[i].quit = 1;
    pthread_barrier_wait(&ThreadBarrierForking);
    for(i=1;i<PARAMS(Threads);i++) pthread_join(ThreadData[i].obj,NULL);
}
#else /* ! USETHREADS */
#define ThreadNo                1 /* number of threads */
#endif /* USETHREADS */

/************************************************************************
*
*    M E M O R Y   M A N A G E M E N T 
*
*************************************************************************
*
* Memory is arranged in slots; a slot consists of several blocks of
* given blocksize. When reallocating, either the number of blocks, the
* blocksize, or both can change. When blocksize changes, the old content
* is adjusted to match the new blocksize.
* Content in the main slots are needed for the duration of the algorithm,
* and the blocksize can change. Temporary slots are used within a single
* iteration, and blocksize is fixed.
*
* void report_memory_usage(channel,prompt)
*   report the cnages since the last report to the given channel
*
* type *get_memory_ptr(type,slot)
*   retrieve the first memory block in the given slot
*
* void clear_memory_slote(void)
*   initialize all memory slots, should be called first
*
* void yalloc(type,slot,blocks,blocksize)
*   initialize a main slot; each block has blocksize words of type
*
* void yrequest(type,slot,blocks,blocksize)
*   request memory allocation change in the given slot. Should be
*   followed by calling reallocmem()
*
* int reallocmem(void)
*   perform previously requested memory allocation changes.
*
* void yfree(slot)
*   release all memory from the given slot
*
* void talloc(type,slot,blocks,blocksize)
*   initialize a temporary slot; previous content is released first.
*   use only for global memory allocation
*
* void talloc2(type,slot,threadID,blocks,blocksize)
*   same as talloc for private memory allocation
*
* void trequest(slot,threadID,blocks)
*   reallocate private memory to 'blocks' many blocks
*/

typedef enum {	/* main memory slots */
M_VertexCoordStore	= 0,	/* vertex coordinates; VertexCoordStore */
M_FacetCoordStore,		/* facet coordinates; FacetCoordStore */
M_VertexAdjStore,		/* adjacency list of vertices; VertexAdjStore*/
M_FacetAdjStore,		/* adjacency list of facets; FacetAdjStore */
M_FacetLiving,			/* facet bitmap; FacetLiving */
M_FacetFinal,			/* facet bitmap; FacetFinal */
M_MAINSLOTS,			/* last main slot index */
		/* temporary slots - threads read these */
M_FacetDistStore=M_MAINSLOTS,	/* facet distances; FacetDistStore */
M_FacetPosnegList,		/* indices of positive/negative facets; FacetPosnegList */
		/* threads write these */
M_THREAD_SLOTS,
M_VertexList=M_THREAD_SLOTS,	/* vertex indices in the intersection of two facets */
M_FacetWork,			/* facet bitmap; FacetWork */
M_VertexArray,			/* calculating facet equations; VertexArray */
M_NewFacetCoordStore,		/* new facet coordinates */
M_NewFacetAdjStore,		/* adjacency list of new facets */

M_THREAD_RESERVED_PLACE,	/* memory block for the next thread */
M_THREAD_RESERVED_PLACE_END = M_THREAD_SLOTS+((M_THREAD_RESERVED_PLACE-M_THREAD_SLOTS)*MAX_THREADS)-1,

M_MSLOTSTOTAL			/* total number of managed memory */
} memslot_t;

/* int NUM_M_THREAD_SLOTS
*    number of thread-specific memory slots */
#define NUM_M_THREAD_SLOTS	(M_THREAD_RESERVED_PLACE-M_THREAD_SLOTS)

/* memslot_t M_thread(slot,threadId)
*    memory slot id assigned to the given thread */
#define M_thread(slot,threadId)		\
    ((slot)+(threadId)*NUM_M_THREAD_SLOTS)

/* slot title for reporting */
#define TM_VertexCoordStore	"VertexCoord"
#define TM_FacetCoordStore	"FacetCoord"
#define TM_VertexAdjStore	"VertexAdj"
#define TM_FacetAdjStore	"FacetAdj"
#define TM_FacetLiving		"FacetLiving"
#define TM_FacetFinal		"FacetFinal"
#define TM_FacetDistStore	"FacetDist"
#define TM_FacetPosnegList	"FacetSign"
#define TM_VertexList		"VertexList"
#define TM_FacetWork		"FacetWork"
#define TM_VertexArray		"VertexArray"
#define TM_NewFacetCoordStore	"NewFacetCoord"
#define TM_NewFacetAdjStore	"NewFacetAdj"

/* MEMSLOT
*    the memory slot structure */

typedef struct { /* memory slot structure */
  size_t blocksize;	/* block size (in bytes) */
  size_t blockno;	/* number of blocks */
  size_t newblocksize;	/* new block size */
  size_t newblockno;	/* new block count */
  size_t rsize;		/* real size in bytes */
  void   *ptr;		/* the actual value */
  const char *title;	/* block title */
  const char *type;	/* basic type */
  size_t bsize;		/* size of item type */
  size_t rreport;	/* rsize at last report */
} MEMSLOT_t;

/* struct MEMSLOT_t memory_slots[]
*    static array containing for each slot the actual blocksize, number
*    of blocks, and a pointer to the actual location. The location can
*    change when reallocating any other memory slot.
*/
static MEMSLOT_t memory_slots[M_MSLOTSTOTAL]; /* memory slots */

/* bool OUT_OF_MEMORY
*    flag indicating whether we are out of memory.
*/
#define OUT_OF_MEMORY	dd_stats.out_of_memory

/* type *get_memory_ptr(type,slot)
*    the actual memory block in the given slot.
*      type:  storage type of the items in the block
*      slot:  the memory slot
*/
#define get_memory_ptr(type,slot)	\
    ((type *)memory_slots[slot].ptr)

/* void report_memory_usage(channel,prompt)
*    for each used slot report the blocksize, number of blocks,
*    total memory used by this slot, and the actual pointer.
*      channel: report type channel
*      force:   report even if there is not change
*      prompt:  prompt for the first line
*/
void report_memory_usage(report_type ch,int force,const char *prompt)
{int i; MEMSLOT_t *ms; char buff[50];
    for(i=0,ms=&memory_slots[0];i<M_MSLOTSTOTAL;i++,ms++) 
    if(ms->ptr && ms->bsize>0 && (force!=0 || ms->rsize != ms->rreport) ) {
        ms->rreport=ms->rsize;
        if(prompt){
            report(ch,"%s\n          ---slot----------+--size---+---blocks--+---blocksize---\n",prompt);
            prompt=NULL;
        }
        if(ms->rsize<1000ul){ sprintf(buff,"%zu",ms->rsize); }
        else if(ms->rsize<1000000ul){ sprintf(buff,"%.2fk",(double)(ms->rsize)*1e-3); }
        else if(ms->rsize<1000000000ul){ sprintf(buff,"%.2fM",(double)(ms->rsize)*1e-6); }
        else { sprintf(buff,"%.2fG",(double)(ms->rsize)*1e-9); }
        report(ch,"          %2d %-14s|%8s |%10zu | %zu * %s\n",
           i+1,ms->title,buff,ms->blockno,ms->blocksize/ms->bsize,ms->type);
    }
    if(prompt==NULL){
       report(ch,"          -----------------+---------+-----------+---------------\n");
    }
}

/* void clear_memory_slots(void)
*    clear all entries in all slots */
#define clear_memory_slots()    \
    memset(&memory_slots[0],0,sizeof(memory_slots))

/* void yalloc(type,slot,n,bsize)
*    initializes a main memory slot by requesting n blocks; each
*    block is an array of bsize elements of the given type.
       type:  storage type of an item
       slot:  memory slot
       n:     number of blocks
       bsize: number of items in a block
*/
#define yalloc(type,slot,n,bsize)	\
    AUX_init_main_slot(slot,n,bsize,sizeof(type),T##slot,mkstringof(type))

/* void AUX_init_main_slot(slot,nno,n,bsize,title,type)
*    initializes a main slot by allocating the requested memory.
*     slot:    memory slot to be initialized
*     nno:     number of initial blocks
*     n:       number of items in a block
*     bsize:   size of an item (in bytes)
*     title:   memory slot title (string) for reporting
*     type:    item type (string) for reporting
*/
static void AUX_init_main_slot(memslot_t slot, size_t nno, size_t n, size_t bsize,
                           const char *title, const char *type)
{size_t total,nsize; MEMSLOT_t *ms;
    if(OUT_OF_MEMORY) return;
    ms=&memory_slots[slot];
    ms->newblocksize=0;
    ms->newblockno=0;
    ms->rreport=0;
    ms->bsize=bsize;	// item size in bytes
    ms->title=title;	// slot title
    ms->type=type;	// item type as string
    nsize=n*bsize;	// block size in bytes
    total=nno*nsize;	// total requested size in bytes
    ms->blocksize=nsize;
    ms->blockno=nno;	// number of requested blocks
    ms->rsize=total;
    ms->ptr=malloc(total);
    if(!ms->ptr){
        report(R_fatal,
           "Out of memory for slot=%d (%s), blocksize=%zu, n=%zu\n",
           slot,title,nsize,nno);
        OUT_OF_MEMORY=1;
        return;
    }
    dd_stats.total_memory += total;
    if(dd_stats.total_memory>dd_stats.max_memory)
             dd_stats.max_memory=dd_stats.total_memory;
    memset(ms->ptr,0,total);
    return;
}

/* void yrequest(type,slot,n,bsize)
*    request memory at a main slot; should be followed by calling
*    reallocmem().
*       type:  storage type of an item
*       slot:  memory slot
*       n:     number of blocks requested
*       bsize: number of items in a block
*/
#define yrequest(type,slot,n,bsize)	\
    AUX_request_main_mem(slot,n,(bsize)*sizeof(type))

/* void AUX_request_main_mem(slot,nno,nsize)
*    records the requested block count and block size. The actual 
*    memory allocation is done by reallocmem()
*      slot:  memory slot
*      nno:   number of blocks requested
*      nsize: size of a block (in bytes)
*/
static inline void AUX_request_main_mem(memslot_t slot, size_t nno, size_t nsize)
{MEMSLOT_t *ms;
    ms=&memory_slots[slot];
    if(ms->blocksize==nsize && nno<=ms->blockno) return;
    ms->newblocksize=nsize;
    ms->newblockno=nno;
}    

/* int reallocmem(void)
*    allocate previously requested memory; if successful, adjust
*    blocks to the given count and size, and clear the new part.
*    Return 1 if out of memory, otherwise return 0. */
static int reallocmem(void)
{MEMSLOT_t *ms; int j,success; size_t total; void *ptr;
    for(j=0,success=1,ms=&memory_slots[0];
        success && j<M_MAINSLOTS; j++,ms++
    ){
        if(ms->newblocksize){
            total=ms->newblocksize*ms->newblockno;
            if(ms->rsize<total){
                dd_stats.memory_allocated_no++;
                ptr=realloc(ms->ptr,total);
                if(ptr){
                    ms->ptr=ptr; 
                    dd_stats.total_memory += total-ms->rsize;
                    if(dd_stats.total_memory>dd_stats.max_memory)
                          dd_stats.max_memory=dd_stats.total_memory;
                    ms->rsize=total;
                }
                else { success=0; }
            }
        }
    }
    if(!success){ OUT_OF_MEMORY=1; return 1; }
    // adjust block structure
    for(j=0,ms=&memory_slots[0]; j<M_MAINSLOTS; j++,ms++){
        if(ms->newblocksize){
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
    // shrink too large allocations
            total=ms->blocksize*ms->blockno;
            if(ms->rsize>total+DD_HIGHWATER){
                ptr=realloc(ms->ptr,total+DD_LOWWATER);
                if(ptr){
                    ms->ptr=ptr;
                    dd_stats.total_memory -= ms->rsize-(total+DD_LOWWATER);
                    ms->rsize=total+DD_LOWWATER;
                }
            }
        }
    }
    return 0;
}

/* void yfree(slot)
*    free the memory in the given slot */
inline static void yfree(memslot_t slot)
{MEMSLOT_t *ms;
    ms=&memory_slots[slot];
    ms->newblocksize=0;
    ms->newblockno=0;
    ms->blockno=0;
    ms->rsize=0;
    if(ms->ptr) free(ms->ptr);
    ms->ptr=(void*)0;
}

/* void talloc(type,slot,n,bsize)
*    request initial memory for a temporary slot. The allocated memory
*    is not cleared
*        type:  storage type of an item
*        slot:  memory slop
*        n:     number of requested blockss
*        bsize: size of an item in bytes
*  void talloc2(type,slotname,slot,n,bsize)
*    talloc() with explicit slot name */
#define talloc(type,slot,n,bsize)	talloc2(type,slot,0,n,bsize)
#define talloc2(type,slot,threadId,n,bsize)	\
    AUX_init_temp_slot(M_thread(slot,threadId),n,bsize,\
         sizeof(type),T##slot,mkstringof(type))

/* void AUX_init_temp_slot(slot,nno,n,bsize,title,type)
*    requests nno blocks, each of size n*bsize at the given slot.
*    The memory is not cleared; should check OUT_OF_MEMORY
*      slot:    memory slot to be initialized
*      nno:     number of requested blocks
*      n:       number of items in each block
*      bsize:   size of an item in bytes
*      title:   memory slot title (string) for reporting
*      type:    item type (string) for reporting 
*/
static void AUX_init_temp_slot(memslot_t slot, size_t nno, size_t n, size_t bsize,
                           const char *title, const char *type)
{size_t total,nsize; MEMSLOT_t *ms;
    if(OUT_OF_MEMORY) return;
    ms=&memory_slots[slot];
    nsize=n*bsize;
    if(ms->bsize!=bsize || ms->blocksize!=nsize || ms->blockno<nno){
        ms->bsize=bsize; ms->blocksize=nsize; ms->blockno=nno;
    }
    ms->title=title; ms->type=type;
    total=nno*nsize;
    if(total <= ms->rsize) return;
    if(ms->ptr){ 
        free(ms->ptr);
        dd_stats.total_memory -= ms->rsize;
        dd_stats.memory_allocated_no++;
    }
    ms->rsize=total;
    dd_stats.total_memory += total;
    if(dd_stats.total_memory>dd_stats.max_memory)
          dd_stats.max_memory=dd_stats.total_memory;
    ms->ptr = malloc(total);
    if(!ms->ptr){
        report(R_fatal,"Out of memory for slot=%d (%s), blocksize=%zu, n=%zu\n",
            slot,title,nsize,nno);
        OUT_OF_MEMORY=1;
    }
}

/* void trequest(slot,n)
*    expand the temporary slot to n blocks, keep the previous block
*    size. The new memory is not cleared. */
#define trequest(slot,threadId,n)		\
    AUX_request_temp_mem(M_thread(slot,threadId),n)

/* void AUX_request_temp_mem(slot,nno)
*    request more blocks for the initialized temporary memory slot
*       slot:  memory slot
*       nno:   number of blocks required
*/
static inline void AUX_request_temp_mem(memslot_t slot,size_t nno)
{size_t total; MEMSLOT_t *ms; void *ptr;
    if(OUT_OF_MEMORY) return;
    ms=&memory_slots[slot];
    total = ms->blocksize*nno;
    if(total<=ms->rsize) return;
    dd_stats.memory_allocated_no++;
    ptr=realloc(ms->ptr,total);
    if(!ptr){ OUT_OF_MEMORY=1; return; }
    dd_stats.total_memory += total - ms->rsize;
    if(dd_stats.total_memory>dd_stats.max_memory)
          dd_stats.max_memory=dd_stats.total_memory;
    ms->rsize=total;
    ms->blockno=nno;
    ms->ptr=ptr;
}

/* void free_adjacency_lists(void)
*    after a break request, vertex and facet adjancy lists are not used
*    release the memory */
void free_adjacency_lists(void)
{   yfree(M_VertexAdjStore); yfree(M_FacetAdjStore); }

/**********************************************************************
*
*     B I T M A P S
*
*************************************************************************
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
*
* bool extract_bit(from,cnt)
*    extract bit cnt from a BITMAP_t array
*
* void clear_bit(from,cnt)
*    clear bit cnt from a BITMAP_t array
*
* void set_bit(from,cnt)
*    set bit cnt in the BITMAP_t array
*
* int get_bitcount(BITMAP_t v)
*    count the bits set in a single BITMAP_t word
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

/* bool extract_bit(BITMAP_t bm[], int cnt)
*  void clear_bit  (BITMAP_t bm[], int cnt)
*  void set_bit    (BITMAP_t bm[], int cnt) */
#define extract_bit(bm,cnt)     \
    (((bm)[(cnt)>>packshift]>>((cnt)&packmask))&1)

#define clear_bit(bm,cnt)       \
    (bm)[(cnt)>>packshift] &= ~(BITMAP1<<((cnt)&packmask))

#define set_bit(bm,cnt)         \
    (bm)[(cnt)>>packshift] |= BITMAP1<<((cnt)&packmask)

/* int get_bitcount(BITMAP_t v)
*    count the the number of bits in v
*/
#ifndef NO_ASM
/* this routine uses the popcnt (population count) machine code
   to get the number of bits in a word. */
static inline int get_bitcount(BITMAP_t v)
{register BITMAP_t res;
    asm ("popcnt %[w], %[t]"
         :[t] "=rm" (res)
         :[w] "rm"  (v));
    return (int)res;
}

#else /* no assembler code */
static int get_bitcount(BITMAP_t v)
{
//    these exceptional cases do not seem to help
//    if(v==BITMAP0){ return 0; }
//    if(v==~BITMAP0){ return 1<<packshift; }
#   ifdef BITMAP_32	/* 32 bit bitmap */
    v=(v&0x55555555u)+((v>>1)&0x55555555u);
    v=(v&0x33333333u)+((v>>2)&0x33333333u);
    v=(v&0x0F0F0F0Fu)+((v>>4)&0x0F0F0F0Fu);
    v=(v&0x00FF00FFu)+((v>>8)&0x00FF00FFu);
    return (int)((v&0x0000FFFF)+(v>>16));
#   else		/* 64 bit bitmap */
    v=(v&0x5555555555555555ul)+((v>>1)&0x5555555555555555ul);
    v=(v&0x3333333333333333ul)+((v>>2)&0x3333333333333333ul);
    v=(v&0x0F0F0F0F0F0F0F0Ful)+((v>>4)&0x0F0F0F0F0F0F0F0Ful);
    v=(v&0x00FF00FF00FF00FFul)+((v>>8)&0x00FF00FF00FF00FFul);
    v=(v&0x0000FFFF0000FFFFul)+((v>>16)&0x0000FFFF0000FFFFul);
    return (int)((v&0xFFFF)+(v>>32));
#   endif
}

#endif /* NO_ASM */

/*************************************************************************
*
*    V E R T I C E S    A N D    F A C E T S
*
**************************************************************************
*
* int DIM
*   problem dimension = number of objectives
* int VertexSize, NextVertex, MaxVertices, VertexBitmapBlockSize
* int FacetSize, NextFacet, MaxFacets, FacetBitmapBlockSize
*   size of a blocks; next free block index; maximal available
*   index; size of the corresponding bitmap block
* int ThisVertex
*   index of facet we are adding to the approximation
* double *VertexCoords(vno), BITMAP_t *VertexAdj(vno)
* double *FacetCoords(fno), BITMAP_t *FacetAdj(fno)
*   the memory block and adjacency block of a vertex and a facet
*
* BITMAP_t *FacetLiving, *FacetFinal
*   bitmaps marking valid and final facets
* double FacetDist[0 .. MaxFacets-1]
*   distance of ThisVertex from living facets
* int *FacetPosnegList
*   positive and negative facet indices 
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
*   o  adjacency lists of vertices and facets are stored as bitmaps
*   o  for each facet its distance from the newly added vertex v; this
*      is computed as
*         distance = f[1]*v[i]+f[2]*v[2]+...+f[DIM]*v{DIM]+f[DIM+1];
*      the polytope is on the >=0 size of each facet.
*   o  two attributes (stored as bitmaps) telling whether a facet is
*      a living one (has not been deleted earlier), and whether this
*      facet is known to be the facet of the final polytope.
*
* All indices start at 0, including vertices, facets, coordinates. */

#define DIM	PARAMS(ProblemObjects)

static int
  VertexSize,		// size of the vertex coordinate block
  NextVertex,		// next free slot for a vertex
  ThisVertex,		// vertex we are working on
  MaxVertices,		// upper bound for vertex numbers
  VertexBitmapBlockSize,// size of the vertex bitmap block
  FacetSize,		// size of the facet equation block
  NextFacet,		// next free slot for a facet
  MaxFacets,		// number of bits in a facet bitmap
  FacetBitmapBlockSize;	// size of the facet bitmap block

/* double *VertexCoords(vno); BITMAP_t *VertexAdj(vno);
*  double *FacetCoords(fno);  BITMAP_t *FacetAdj(fno);
*     coordinate and adjacency list of vertex vno and facet fno */
#define VertexCoords(vno)\
    (get_memory_ptr(double,M_VertexCoordStore)+((vno)*VertexSize))
#define VertexAdj(vno)          \
    (get_memory_ptr(BITMAP_t,M_VertexAdjStore)+((vno)*FacetBitmapBlockSize))
#define FacetCoords(fno)        \
    (get_memory_ptr(double,M_FacetCoordStore)+((fno)*FacetSize))
#define FacetAdj(fno)           \
    (get_memory_ptr(BITMAP_t,M_FacetAdjStore)+((fno)*VertexBitmapBlockSize))

/* BITMAP_f *FacetLiving, *FacetFinal; double *FacetDist; int* FacetPosnegList; */
#define FacetLiving		\
    get_memory_ptr(BITMAP_t,M_FacetLiving)
#define FacetFinal		\
    get_memory_ptr(BITMAP_t,M_FacetFinal)
#define FacetDist(fno)  	\
    get_memory_ptr(double,M_FacetDistStore)[fno]
#define FacetPosnegList		\
    get_memory_ptr(int,M_FacetPosnegList)

/* BOOL is_livingFacet(fno)
*     check of bit 'fno' is set in bitmap FacetLiving
* BOOL is_finalFacet(fno)
*     check if bit 'fno' is set in bitmap FacetFinal */

#define is_livingFacet(fno)     	extract_bit(FacetLiving,fno)
#define is_finalFacet(fno)		extract_bit(FacetFinal,fno)

/* void mark_facet_as_final(fno)
*     exported version of set_bit(FacetLiving,fno)
* void clear_facet_as_living(fno)
*     exported version of clear_bit(FacetLiving,fno)
*  void intersect_vertexAdj_FacetLiving(vno)
*     clear bits in the adjacency list which are not living 
*  void copy_facet_to(coords,adj_bitmap,to)
*     move the facet with coords and adjacency map to the given address
*  void clear_VertexAdj(fno); void clear_FacetAdj(vno)
*     clear the adjacency list of facet 'fno' and vertex 'vno'
* void copy_FacetLiving_to(where)
*     copy the FacetLiving bitmap to the given location */

void mark_facet_as_final(int fno)
{   set_bit(FacetFinal,fno); }
void clear_facet_as_living(int fno)
{   clear_bit(FacetLiving,fno); }
inline static void intersect_VertexAdj_FacetLiving(int vno)
{int i;
    for(i=0;i<FacetBitmapBlockSize;i++) VertexAdj(vno)[i] &= FacetLiving[i];
}
inline static void copy_facet_to(double *coords,BITMAP_t *adj,int fno)
{   memcpy(FacetCoords(fno),coords,FacetSize*sizeof(double));
    memcpy(FacetAdj(fno),adj,VertexBitmapBlockSize*sizeof(BITMAP_t)); }
inline static void clear_VertexAdj(int vno)
{   memset(VertexAdj(vno),0,FacetBitmapBlockSize*sizeof(BITMAP_t)); }
inline static void clear_FacetAdj(int fno)
{   memset(FacetAdj(fno),0,VertexBitmapBlockSize*sizeof(BITMAP_t)); }
inline static void copy_FacetLiving_to(BITMAP_t *where)
{   memcpy(where,FacetLiving,FacetBitmapBlockSize*sizeof(BITMAP_t)); }


/*************************************************************************
*
*    R E P O R T I N G    
*
**************************************************************************/

DD_STATS dd_stats;

/* Statistics:
*  void get_dd_facetno(void)
*     compute the number of living and final facets int dd_stats
*  int get_vertexnum(void), get_facetnum(void)
*     returnt the number of added vertices and the number of living facets
*
* Rounding proceudres:
*  double round(double x)
*    if x is closer than ScaleEps to an integer, return that integer
*  bool closetoint(double x)
*    true if x is close thatn ScaleEps to an integer
*  int gcd(int a,int b)
*    the greatest common divisor of non-negative a and b
*  int lcm(a,b)
*    the least common divisor of positive a and b
*  int denum(double x)
*    if d*x is close to an integer for a small integer d, return d
*
* Reporting, retrieving:
*  void print_vertex(report_type channel, double coords[0:DIM-1](
*    report the vertex using fractional or floating point format.
*  void print_vertices(report_type channel)
*    report all stored vertices using print_vertex().
*  void print_facet(report_type channel, int fno)
*    prints the coordinates of the facet using scaling and rounding 
*  void print_facets(report_type channel)
*    report all facets using scaling and rounding format.
*  void get_facet_into(fno, double v[0:dim])
*    store the facet equation to the provided space. The numbers are 
*    scaled and rounded to the closest integer 
*/

void get_dd_facetno(void) // get the number living and final facets
{int i,fno; BITMAP_t vl,vf;
    dd_stats.living_facets_no=-1;
    dd_stats.final_facets_no=-1;
    for(i=0;i<FacetBitmapBlockSize;i++){
        vl=FacetLiving[i];
        dd_stats.living_facets_no+=get_bitcount(vl);
        vf=FacetFinal[i];
        dd_stats.final_facets_no+=get_bitcount(vf);
        if(~vl & vf){ /* consistency checking */
             fno=i<<packshift;
             while(!((~vl & vf)&1)){ fno++; vl>>=1; vf>>=1; }
             report(R_err,"Consistency error: final but not living facet %d\n",fno);
        }
    }
    if(dd_stats.final_facets_no<0) dd_stats.final_facets_no=0;
    if(dd_stats.living_facets_no<0) dd_stats.living_facets_no=0;
}

int get_vertexnum(void)
{   return NextVertex-DIM; }

int get_facetnum(void)
{int i,total; BITMAP_t v;
    total=-1;
    for(i=0;i<FacetBitmapBlockSize;i++){
        v=FacetLiving[i]; total += get_bitcount(v);
    }
    return total;
}

/* Rounding procedures */

#include "round.h" /* round_top() */

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
     round_to(&x);
     return x;
}
inline static int closetoint(double x)
{   if(x<0.0) x=-x;
    x -= (int)(x+0.5); 
    return (x<SCALE_EPS && x>-SCALE_EPS);
}
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
static int denum(double x)
{int i;
    for(i=1;i<1000;i++)if(closetoint(i*x)) return i;
    return 1;
}

/* Reporting, retrieving */

static char *formatvalue(double v)
{static char buf[80]; int d;
    round_to(&v);
    if(!closetoint(v) && PARAMS(VertexAsFraction)){
        d=denum(v); if(closetoint(d*v)){
            sprintf(buf,"%d/%d",(int)(round(d*v)),d);
            return buf;
        }
    }
    sprintf(buf,"%.14g",round(v));
    return buf;
}
void print_vertex(report_type channel, const double v[/*0:DIM-1*/])
{int j; double dir;
    dir = PARAMS(Direction) ? -1.0 : +1.0;
    for(j=0;j<DIM;j++)
        report(channel," %s",formatvalue(dir*v[j]));
    report(channel,"\n");
}
void print_vertices(report_type channel)
{int i;
    for(i=DIM;i<NextVertex;i++){
        report(channel,"V ");
        print_vertex(channel,VertexCoords(i));
    }
}
void print_facet(report_type channel, int fno)
{int j,d; double dir;
    report(channel,"%c ",is_finalFacet(fno)?'F':'f');
    d=1; for(j=0;d<130000 &&j<=DIM;j++){
        d=lcm(d,denum(FacetCoords(fno)[j]));
    }
    dir = PARAMS(Direction) ? -1.0 : +1.0;
    for(j=0;j<DIM;j++)
       report(channel," %.14lg",round(d*FacetCoords(fno)[j]));
    report(channel," %.14lg\n",dir*round(d*FacetCoords(fno)[DIM]));
}
void print_facets(report_type channel)
{int i;
    for(i=1;i<NextFacet;i++)if(is_livingFacet(i))
        print_facet(channel,i);
}

static void get_facet_into(int fno, double *v)
{int j,d;
    d=1; for(j=0;d<130000 && j<=DIM;j++)d=lcm(d,denum(FacetCoords(fno)[j]));
    for(j=0;j<=DIM;j++)v[j]=round(d*FacetCoords(fno)[j]);
}

/*************************************************************************
*
*    I T E R A T I O N
*
**************************************************************************
*
* The main iteration adds a new vertex to the actual approximation. First,
* for each facet 'fno' the distance from the new vertex is stored in
*    double   FacetDist(fno).
* As the next step,
*    int      *FacetPosnegList[0 .. MaxFacets-1] 
* is populated by facet numbers: those with positive distance go the the
* front, those with negative distance go to the end. Then for each 
* positive / negative pair (fp,fn) the following are computed by threads:
*    int      *VertexList(thId)
* vertex numbers adjacent to both (fp,fn);
*    BITMAP_t *FacetWork(thId)
* bitmap of facets adjacent to vertices in VertexList. If (fp,fn) is a
* ridge, it gives a new facet together with ThisVertex
*    double   *NewFacetCoords(thId,fno)
*    BITMAP_t *NewFacetAdj(thId,fno)
* coordinates and the adjacency list of the facet created by thread thId.
* When there are no threads, then thId=0. */

#define VertexList(thId)		\
    get_memory_ptr(int,M_thread(M_VertexList,thId))
#define FacetWork(thId)			\
    get_memory_ptr(BITMAP_t,M_thread(M_FacetWork,thId))
#define NewFacetCoords(thId,fno)	\
    (get_memory_ptr(double,M_thread(M_NewFacetCoordStore,thId))+\
          ((fno)*FacetSize))
#define NewFacetAdj(thId,fno)		\
    (get_memory_ptr(BITMAP_t,M_thread(M_NewFacetAdjStore,thId))+\
          ((fno)*VertexBitmapBlockSize))
#define VertexArray(thId)		\
    get_memory_ptr(double,M_thread(M_VertexArray,thId))

/* data set byt the threads */
static int
    NewFacet[MAX_THREADS],	// number of newly created facets
    MaxNewFacets[MAX_THREADS],	// available space
    ErrorNo[MAX_THREADS];	// consistency errors

/* void move_NewFacet_th(threadId,fno)
*     move NewFacet[threadId]-th element to the slot 'fno' */

inline static void move_NewFacet_th(int thId,int fno)
{   copy_facet_to(NewFacetCoords(thId,NewFacet[thId]),
                  NewFacetAdj(thId,NewFacet[thId]),fno); }

/***********************************************************************
* initialize all data structures and add the very first vertex
*
* int init_dd_structure(int vertexno,int facetno)
*   initialize the DD algorithm structures. Return value:
*  0  initialization is OK
*  1  problem: dimension is too large, or out of memory
* void init_dd(double first_vertex[0:DIM-1]
*   create the initial inner approximation with the given vertex and
*   all ideal vertices; add the ideal facet and DIM many other facets */

int init_dd_structure(int vertexno, int facetno)
{int i,j;
    if(DIM<1||DIM>MAXIMAL_ALLOWED_DIMENSION){
        report(R_fatal,"init_dd: dimension %d is out of allowed range (1..%d)\n",
           DIM,MAXIMAL_ALLOWED_DIMENSION);
           return 1;
    }
    MaxVertices=vertexno;
    if(MaxVertices<DD_INITIAL_VERTEXNO) MaxVertices=DD_INITIAL_VERTEXNO;
    MaxFacets=facetno;
    if(MaxFacets<DD_INITIAL_FACETNO) MaxFacets=DD_INITIAL_FACETNO;
    while(MaxVertices+5<DIM) MaxVertices += (DD_VERTEX_ADDBLOCK<<packshift);
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
    dd_stats.vertexno=0; // the ideal vertex is not counted
    // clear memory slots
    clear_memory_slots();
    yalloc(double,M_VertexCoordStore,MaxVertices,VertexSize); // VertexCoordStore
    yalloc(double,M_FacetCoordStore,MaxFacets,FacetSize); // FacetCoordStore
    yalloc(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize); // VertexAdjStore
    yalloc(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize); // FacetAdjStore
    yalloc(BITMAP_t,M_FacetLiving,1,FacetBitmapBlockSize); // FacetLiving
    yalloc(BITMAP_t,M_FacetFinal,1,FacetBitmapBlockSize); // FacetFinal
    if(OUT_OF_MEMORY) return 1;
    dd_stats.memory_allocated_no=1;
    NextVertex=0; NextFacet=0;
    // the ideal vertices with index 0..DIM-1
    for(i=0;i<DIM;i++){ // the i-th ideal vertex
        for(j=0;j<DIM;j++) VertexCoords(NextVertex)[j]= i==j ? 1.0: 0.0;
        clear_VertexAdj(NextVertex);
        set_bit(VertexAdj(NextVertex),0); // it is on the ideal facet
        NextVertex++;
    }
    // the ideal facet
    for(j=0;j<=DIM;j++) FacetCoords(NextFacet)[j]= j==DIM? 1.0 : 0.0;
    clear_FacetAdj(NextFacet);
    for(j=0;j<DIM;j++)set_bit(FacetAdj(NextFacet),j);
    set_bit(FacetLiving,NextFacet); set_bit(FacetFinal,NextFacet);
    NextFacet++;
    return 0;
}

/* void init_dd(double coords[0..DIM-1])
*    set up the first non-ideal vertex and DIM non-ideal facets */
void init_dd(const double *coords)
{int i,j;
    // set up the first vertex
    for(j=0;j<DIM;j++) VertexCoords(NextVertex)[j]=coords[j];
    clear_VertexAdj(NextVertex);
    // and DIM facets
    for(i=0;i<DIM;i++){
        for(j=0;j<DIM;j++) FacetCoords(NextFacet)[j]= i==j ? 1.0 : 0.0;
        FacetCoords(NextFacet)[DIM]=-coords[i];
        clear_FacetAdj(NextFacet);
        set_bit(FacetAdj(NextFacet),NextVertex);
        set_bit(VertexAdj(NextVertex),NextFacet);
        for(j=0;j<DIM;j++) if(i!=j){
            set_bit(FacetAdj(NextFacet),j); set_bit(VertexAdj(j),NextFacet);
        }
        set_bit(FacetLiving,NextFacet);
        NextFacet++;
    }
    NextVertex++;
    dd_stats.vertexno++;
}

/* int add_initiali_vertex(double coords[0:dim-1])
*    collect all (non-ideal) vertices of the inner approximation first
*  int add_initial_facet(int final,double coords[0:dim])
*    create a consistent approximating polytope
*    return value: 0: OK, 1: some error */
int add_initial_vertex(const double coords[/* 0..DIM-1 */])
{int j; double dir;
    if(NextVertex>=MaxVertices){
        report(R_fatal,"Resume: more vertices than specified (%d)\n",MaxVertices);
        return 1; 
    }
    if(NextFacet>1){
        report(R_fatal,"Resume: vertices should come before facets\n");
        return 1;
    }
    dir = PARAMS(Direction) ? -1.0 : +1.0; 
    for(j=0;j<DIM;j++) VertexCoords(NextVertex)[j]=dir*coords[j];
    clear_VertexAdj(NextVertex);
    NextVertex++;
    dd_stats.vertexno++;
    return 0;
}

int add_initial_facet(int final,const double coords[/* 0..DIM */])
{int i,j; double w;
    if(NextFacet>=MaxFacets){
        report(R_fatal,"Resume: more facets than specified (%d)\n",MaxFacets);
        return 1;
    }
    w=0.0; for(j=0;j<DIM;j++){
        w+=coords[j];
        if(coords[j]<0.0){
            report(R_fatal,"Resume: facet %d has negative coefficient\n",NextFacet);
            return 1;
        }
    }
    if(w<PARAMS(PolytopeEps)){
        report(R_fatal,"Resume: facet %d has all zero coefficients\n",NextFacet);
        return 1;
    }
    w=1.0/w;
    for(j=0;j<=DIM;j++){
        FacetCoords(NextFacet)[j]=w*coords[j];
    }
    set_bit(FacetLiving,NextFacet);
    if(final) set_bit(FacetFinal,NextFacet);
    clear_FacetAdj(NextFacet);
    for(i=0;i<DIM;i++){ // ideal vertices
        if(FacetCoords(NextFacet)[i]<PARAMS(PolytopeEps)){
            set_bit(FacetAdj(NextFacet),i); set_bit(VertexAdj(i),NextFacet);
        }
    }
    for(i=DIM;i<NextVertex;i++){ // other vertices
        w=FacetCoords(NextFacet)[DIM];
        for(j=0;j<DIM;j++){ w += FacetCoords(NextFacet)[j]*VertexCoords(i)[j]; }
        if(w<-PARAMS(PolytopeEps)){
            report(R_fatal,"Resume: vertex %d is on the negative side of facet %d (%lg)\n",i,NextFacet,w);
            return 1;
        }
        if(w<PARAMS(PolytopeEps)){ // adjacent
            set_bit(FacetAdj(NextFacet),i); set_bit(VertexAdj(i),NextFacet);
        }
    }
    NextFacet++;
    return 0;
}

/***********************************************************************
* Request memory to accomodate more vertices and facets
*
* void allocate_vertex_block(int bitmap)
*   add addition DD_VERTEX_ADDBLOCK<<packshift more vertices. If bitmap
*   is NOT set, don't allocate bitmaps
* void allocate_facet_block(count)
*   add space for 'count' more facets */

static void allocate_vertex_block(int bitmap)
{   // extend Vertices and FacetBitmap to accommodate more stuff
    if(OUT_OF_MEMORY) return;
    dd_stats.vertices_allocated_no ++;
    dd_stats.vertices_allocated += (DD_VERTEX_ADDBLOCK<<packshift);
    MaxVertices += (DD_VERTEX_ADDBLOCK<<packshift);
    // tell the memory handling part how much storage space we would need.
    yrequest(double,M_VertexCoordStore,MaxVertices,VertexSize);
    if(bitmap){
       VertexBitmapBlockSize += DD_VERTEX_ADDBLOCK;
       yrequest(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize);
       yrequest(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize);
    }
    // and do the reallocation
    if(reallocmem()){  // out of memory, don't increase the values
        MaxVertices -= (DD_VERTEX_ADDBLOCK<<packshift);
        if(bitmap) VertexBitmapBlockSize -= DD_VERTEX_ADDBLOCK;
    }
}

static void allocate_facet_block(int count)
{int total;
    if(OUT_OF_MEMORY) return;
    for(total=DD_FACET_ADDBLOCK<<packshift;total<count; 
            total += DD_FACET_ADDBLOCK<<packshift);
    dd_stats.facets_allocated_no ++;
    dd_stats.facets_allocated += total;
    MaxFacets += total;
    FacetBitmapBlockSize = (MaxFacets+packmask)>>packshift;
    yrequest(double,M_FacetCoordStore,MaxFacets,FacetSize);
    yrequest(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_FacetLiving,1,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_FacetFinal,1,FacetBitmapBlockSize);
    if(reallocmem()){ // out of memory
        MaxFacets -= total;
        FacetBitmapBlockSize = (MaxFacets+packmask)>>packshift;
    }
}

/***********************************************************************
* int facet_intersection(f1,f2)
*    intersect the vertex adjacency lists of f1 and f2; return the
*    number of vertices adjacent to both f1 and f2 */
inline static int facet_intersection(int f1, int f2)
{int i,total; register BITMAP_t *L1,*L2;
    total=0; L1=FacetAdj(f1); L2=FacetAdj(f2);
    for(i=0;i<VertexBitmapBlockSize;i++,L1++,L2++)
        total += get_bitcount((*L1)&(*L2));
    return total;
}

/***********************************************************************
* int get_new_facetno(int threadId)
*    return the next available facet number, asking for memory if necessary.
*    The return value is -1 if cannot allocate more memory.
*/

inline static int get_new_facetno(int threadId)
{int i;
    i=NewFacet[threadId];
    if(i>=MaxNewFacets[threadId]){ // no more space, ask memory
        MaxNewFacets[threadId] += DD_FACET_ADDBLOCK<<packshift;
        trequest(M_NewFacetCoordStore,threadId,MaxNewFacets[threadId]);
        trequest(M_NewFacetAdjStore,threadId,MaxNewFacets[threadId]);
        if(OUT_OF_MEMORY) {
            MaxNewFacets[threadId] -= DD_FACET_ADDBLOCK<<packshift;
            return -1;
        }
    }
    NewFacet[threadId]++; // we have one more
    return i;
}

/* double vertex_distance(double *coords, int fno)
*    compute the distance of the vertex given as the first argument
*    from the facet given with its number */
static double vertex_distance(const double *coords, int fno)
{double d=0.0; int i; double *fcoords;
    fcoords=FacetCoords(fno);
    for(i=0;i<DIM;i++){
       d+= (*coords)*(*fcoords);
       coords++; fcoords++;
    }
    d += (*fcoords);
    return d;
}

/* int probe_vertex(double *coords)
*    return the number of negative facets. */
int probe_vertex(double *coords)
{int fno,negfacets;
    dd_stats.probevertex++; negfacets=0;
    for(fno=0;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        if(vertex_distance(coords,fno) < -PARAMS(PolytopeEps)) negfacets++;
    }
    return negfacets;
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
        report(R_err,"solve_lineq: rank is too large, increase LineqEps=%lg\n",
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

/* void recalculate_facet(int fno, BITMAP_t *adj, double *coords, memslot_t tmp)
*    calculate facet equation form the vertices adjacent to it. Complain
*    if out of memory, the system is degenerate, or the old and new coeffs
*    are too far from each other. */
static void recalculate_facet(int fno,BITMAP_t *adj,double *old,int threadId)
{double v,s; int i,j,vno; int an,DIM1; BITMAP_t fc; double *VA;    
    // collect all vertices adjacent to facet fno
    DIM1=DIM+1;
    an=0; for(i=0;i<VertexBitmapBlockSize;i++) an+=get_bitcount(adj[i]);
    // request memory
    talloc2(double,M_VertexArray,threadId,an,DIM+1);
    if(OUT_OF_MEMORY) return;
    VA=VertexArray(threadId); an=0; vno=0;
#define A(i,j)	VA[(i)*DIM1+(j)]
    for(i=0;i<VertexBitmapBlockSize;i++){
        j=vno; fc=adj[i]; while(fc){
            while((fc&7)==0){fc>>=3; j+=3; }
            if(fc&1){ /* store vertex j to A */
                memcpy(&(A(an,0)),VertexCoords(j),DIM*sizeof(double));
                A(an,DIM)= j<DIM ? 0.0 : 1.0 ;
                an++;
            }
            j++; fc>>=1;
        }
        vno += (1<<packshift);
    }
    // all vertices are collected, compute facet equation
    if(an<DIM){
        report(R_err,"recalculate: facet %d has only %d adjacent vertices\n",fno,an);
        dd_stats.numerical_error++;
        return;
    }
    if(solve_lineq(an,VA)){
        report(R_err,"recalculate: adjacency vertex list of facet %d is degenerate\n",fno);
        dd_stats.numerical_error++;
        return;
    }
    /* result is in VA[0..DIM] */
    s=0.0; for(i=0;i<DIM;i++) s+= VA[i]; s=1.0/s;
    for(i=0;i<=DIM;i++){
        VA[i]*=s;
        v=old[i]-VA[i];
        if(v>PARAMS(FacetRecalcEps) || v < -PARAMS(FacetRecalcEps)){
            report(R_warn,"recalculate: numerical instability at facet %d, coord %d (%lg)\n", fno,i,v);
            dd_stats.instability_warning++;
        }
        old[i]=VA[i];
    }
#undef A
}

/* void recalculate_facets(void)
*    go over all facets, and recalculate their equations. 
*  voiod thread_recaluclate(threadId)
*    split all cases int ThreadNo pieces; each thread excutes one of
*    them. If there are no threads then ThreadNo=1 */

static void thread_recalculate(int threadId) // Id goes from 0 to ThreadNo-1
{int fno,step;
    step=ThreadNo;
    for(fno=1+threadId;fno<NextFacet;fno+=step) if(is_livingFacet(fno))
        recalculate_facet(fno,FacetAdj(fno),FacetCoords(fno),threadId);
}

void recalculate_facets(void)
#ifdef USETHREADS
{ thread_execute(thread_recalculate); }
#else /* ! USETHREADS */
{ thread_recalculate(0); }
#endif /* USETHREADS */

/***********************************************************************
* The combinatorial Double Description method
*
* int is_ridge(f1,f2,threadId)
*    combinatorial test to check whether the intersection of f1 and
*    f2 forms a ridge. If there are <DIM-1 vertices both in f1 and
*    f2 then it is not a ridge. Otherwise store indices adjacent to both
*    f1 and f2 in VertexList, then take the intersection of living facets
*    and the facet adjacency lists in VertexList. (f1,f2) is a ridge if
*    and only if the intersection contains f1 and f2 only.
* void create_new_facet(f1,f2,threadId)
*    create a new facet through the ridge formed by f1 and f2, and
*    the vertex ThisVertex. f1*v is negative, f2*v is positive, and these
*    values are stored in FacetDist[]. Recalculate the equation when
*    ExactFacetEq parameter is set 
 void search_ridges(void)
*    go over all facet pairs (f1,f2), f1<0, f2>0 and check if it is a ridge.
*    There are dd_stats.facet_neg negative and dd_stats.facet_pos positive
*    facets stored at FacetPosnegList (from the end, from the beginning)
*    If ridge, call create_new_facet().
* void thread_search_ridges(threadId)
*    split all cases into ThreadNo pieces, each thread executes one of them.
*    It there are no threads, ThreadNo=1, and threadId=0 
*/

inline static int is_ridge(int f1,int f2,int threadId)
{int vertexno,i,j,vlistlen; BITMAP_t v,*va0,*va1;
    if(facet_intersection(f1,f2) < DIM-1)
         return 0; // no - happens oftern
    /* all living facets except for f1 and f2 */
    copy_FacetLiving_to(FacetWork(threadId));
    clear_bit(FacetWork(threadId),f1); clear_bit(FacetWork(threadId),f2);
    /* store vertex indices adjacent to both f1 and f2 in VertexList */
    vertexno=0; vlistlen=0;
    for(i=0;i<VertexBitmapBlockSize;i++){
        if((v=FacetAdj(f1)[i] & FacetAdj(f2)[i])){
            j=vertexno;
            while(v){
                while((v&7)==0){ j+=3; v>>=3; }
                if(v&1){ VertexList(threadId)[vlistlen]=j; vlistlen++; }
                j++; v>>=1;
            }
        }
        vertexno += (1<<packshift);
    }
    /* now we have all vertices adjacent to f1 and f2 in VertexList
       check if the intersection of the adjacency lists of these vertices
       contain only f1 and f2 and nothing else */
    va0=VertexAdj(VertexList(threadId)[0]);
    va1=VertexAdj(VertexList(threadId)[1]);
    for(i=0;i<FacetBitmapBlockSize;i++) if(
       (v=FacetWork(threadId)[i] & va0[i]) && (v &=va1[i])){
          for(j=2;j<vlistlen;j++) v &= VertexAdj(VertexList(threadId)[j])[i];
          if(v) return 0; // no
    }
    return 1; // yes
}

/* void create_new_facet(f1,f2,threadId)*/
inline static void create_new_facet(int f1, int f2, int threadId)
{int newf; double d1,d2,d; int i;
    newf=get_new_facetno(threadId);
    if(newf<0) return; // no memory
    // the adjacency list is the intersection of that of f1 and f2 plus the new vertex
    for(i=0;i<VertexBitmapBlockSize;i++)
        NewFacetAdj(threadId,newf)[i] = FacetAdj(f1)[i] & FacetAdj(f2)[i];
    set_bit(NewFacetAdj(threadId,newf),ThisVertex);
    // compute the coefficients, f1<0, f2 >0
    if(f2==0){ // ideal facet
        for(i=0;i<=DIM;i++){
            NewFacetCoords(threadId,newf)[i]=FacetCoords(f1)[i];
        }
        NewFacetCoords(threadId,newf)[DIM] -= FacetDist(f1);
    } else {   // f1(v)=-d1, f2(v)=d2; (d2*f1 + d1*f2)(v)=0.0
        d1=-FacetDist(f1); d2=FacetDist(f2);
        d=1.0/(d1+d2); d1 *= d; d2 *= d;
        for(i=0;i<=DIM;i++){
            NewFacetCoords(threadId,newf)[i]=
              d2*FacetCoords(f1)[i]+d1*FacetCoords(f2)[i];
        }
    }
    if(PARAMS(ExactFacetEq))
        recalculate_facet(MaxFacets+newf,  // report this number if error
            NewFacetAdj(threadId,newf),    // adjacency list
            NewFacetCoords(threadId,newf), // old coordinates, replaced
            threadId);                     // thread
}

/* void thread_search_ridges(threadId) */
static void thread_search_ridges(int threadId) // threadId goes from 0 to ThreadNo-1
{int i,j,f1,step; int *PosIdx, *NegIdx;
    step=ThreadNo;
    NegIdx = FacetPosnegList+(MaxFacets-1-threadId);
    for(j=threadId;j<dd_stats.facet_neg;j+=step,NegIdx -= step){
       f1=*NegIdx; PosIdx=FacetPosnegList;
       for(i=0;i<dd_stats.facet_pos;i++,PosIdx++){
           if(is_ridge(f1,*PosIdx,threadId))
               create_new_facet(f1,*PosIdx,threadId);
       }
    }
}

/* void search_ridges() */
static void search_ridges(void)
#ifdef USETHREADS
{   thread_execute(thread_search_ridges); }
#else /* !USETHREADS */
{   thread_search_ridges(0); }
#endif /* USETHREADS */


/***********************************************************************
* void make_facet_living(fno)
*    set the "living" flag for this facet, and add it to the adjacency
*    list of all vertices it is adjacent to.
* void compress_from(fno)
*    move living facets into the unoccupied slots starting from 'fno'
* void compress_facets(void)
*    if facet numbers decreased, decrease the storage as well
* void reaquest_main_loop_memory(threadId)
*    allocate temporary memory used by a thread in the main loop
*    VertexList:      list of vertices adjacent to facets f1 and f2
*    FacetWork:       bitmap of facets adjacent to all vertices in VertexList
*    NewFacetCoord:   coordinates of new facet
*    NewFacetAdj:     adjacency bitmap of the new facet 
* void add_new_vertex(double vertex[0:dim-1])
*    add a new vertex which is outside the convex hull of the present
*    approximation. Split facets into positive, negative and zero sets
*    depending where the new vertex is. For each pair of positive/
*    negative facets, check if it is a ridge; if yes, add a new facet.
*    Then throw away negative facets, and add the newly created facets
*    to the approximation.
*/

/* void make_facet_living(facetno) */
static void make_facet_living(int fno)
{int vno,i,j; BITMAP_t fc;
    set_bit(FacetLiving,fno);
    vno=0; for(i=0;i<VertexBitmapBlockSize;i++){
        j=vno; fc=FacetAdj(fno)[i];
        while(fc){
            while((fc&7)==0){ j+=3; fc>>=3; }
            if(fc&1){set_bit(VertexAdj(j),fno);}
            j++; fc>>=1;
        }
        vno+=(1<<packshift);
    }
}

/* void compress_from(int fno) */
static void compress_from(int fno)
{int i;
 fill_holes:
    while(fno<NextFacet && is_livingFacet(fno) ) fno++;
    while(fno<NextFacet  && !is_livingFacet(NextFacet-1) ) NextFacet--;
    if(fno < NextFacet){ // fno is empty, NextFacet-1 is used
        NextFacet--;
        copy_facet_to(FacetCoords(NextFacet),FacetAdj(NextFacet),fno);
        make_facet_living(fno);
        clear_bit(FacetLiving,NextFacet);
        if(is_finalFacet(NextFacet)) set_bit(FacetFinal,fno);
        fno++; goto fill_holes;
    }
    // clear FacetFinal
    for(i=0;i<FacetBitmapBlockSize;i++) FacetFinal[i] &= FacetLiving[i];
    // clear the adjacency list of vertices
    for(i=0;i<NextVertex;i++) intersect_VertexAdj_FacetLiving(i);
}

/* void compress_facets(void) */
static void compress_facets(void)
{int i,fno;
    if(MaxFacets <= get_facetnum()+2*(DD_FACET_ADDBLOCK<<packshift)) return;
    dd_stats.facets_compressed_no++;
    // clear the adjacency list of vertices
    for(i=0;i<NextVertex;i++) intersect_VertexAdj_FacetLiving(i);
    // find the first free facet
    for(i=fno=0;fno<NextFacet && (~FacetLiving[i])==0;i++,fno+=(1<<packshift));
    compress_from(fno);
    // compute new limits
    while(MaxFacets > NextFacet+2*(DD_FACET_ADDBLOCK<<packshift)){
        MaxFacets -= DD_FACET_ADDBLOCK<<packshift;
    }
    FacetBitmapBlockSize = (MaxFacets+packmask)>>packshift;
    yrequest(double,M_FacetCoordStore,MaxFacets,FacetSize);
    yrequest(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_FacetLiving,1,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_FacetFinal,1,FacetBitmapBlockSize);
    if(reallocmem()){ // fatal error
        report(R_fatal,"Compress facets error, problem aborted\n");
        exit(4);      // unrecoverable error
    }
}

/* void request main loop memory(threadId) */
inline static void request_main_loop_memory(int threadId)
{   talloc2(int,M_VertexList,threadId,MaxVertices,1);
    talloc2(BITMAP_t,M_FacetWork,threadId,1,FacetBitmapBlockSize);
    talloc2(double,M_NewFacetCoordStore,threadId,DD_INITIAL_FACETNO,FacetSize);
    talloc2(BITMAP_t,M_NewFacetAdjStore,threadId,DD_INITIAL_FACETNO,VertexBitmapBlockSize);
    NewFacet[threadId]=0;
    MaxNewFacets[threadId]=DD_INITIAL_FACETNO;
}

/* add a new vertex to the approximation */
void add_new_vertex(const double *coords)
{int i,j,fno,threadId,AllNewFacets; BITMAP_t fc; double d;
 int *PosIdx, *NegIdx;
    dd_stats.iterations++; dd_stats.vertexno++;
    if(NextVertex >= MaxVertices){
        compress_facets();
        allocate_vertex_block(1);
    }
    // memory blocks for the outer loop
    talloc(double,M_FacetDistStore,MaxFacets,1);
    talloc(int,M_FacetPosnegList,MaxFacets,1);
    for(threadId=0;threadId<ThreadNo;threadId++)
        request_main_loop_memory(threadId);
    if(OUT_OF_MEMORY){
        dd_stats.data_is_consistent=1;
        return;
    }
    ThisVertex=NextVertex; NextVertex++;
    for(i=0;i<DIM;i++) VertexCoords(ThisVertex)[i]=coords[i];
    clear_VertexAdj(ThisVertex); // clear the adjacency list
    dd_stats.facet_pos=0; dd_stats.facet_neg=0; dd_stats.facet_zero=0;
    PosIdx = FacetPosnegList; // this goes ahead
    NegIdx = FacetPosnegList+MaxFacets; // this goes backward
    for(fno=0;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        d=FacetDist(fno)=vertex_distance(coords,fno);
        if(d>PARAMS(PolytopeEps)){ // positive side 
            *PosIdx = fno; ++PosIdx;
            dd_stats.facet_pos++;
        } else if(d<-PARAMS(PolytopeEps)){ // negative side
            if(is_finalFacet(fno)){
                report(R_warn,"Final facet %d is on the negative side of "
                    "vertex %d (d=%lg)\n",fno,ThisVertex-DIM+1,d);
                dd_stats.instability_warning++;
// it seems to be better to revoke the "final" flag from this facet
                clear_bit(FacetFinal,fno);
            }
            --NegIdx; *NegIdx=fno;
            dd_stats.facet_neg++;
        } else { // this is adjacent to our new vertex
            set_bit(VertexAdj(ThisVertex),fno);
            set_bit(FacetAdj(fno),ThisVertex);
            dd_stats.facet_zero++;
        }
    }
    if(dd_stats.facet_neg==0){ // this vertex is inside the polytope
        dd_stats.facet_new=0;  // no new facet is added
        if(dd_stats.facet_zero<DIM){
            report(R_err,"Next vertex is inside the approximation\n");
            dd_stats.numerical_error++;
        } else { // check whether this is a duplicate vertex
            dd_stats.instability_warning++;
            for(j=ThisVertex-1; j>=DIM && j>ThisVertex-100; j--){
              for(i=0;i<DIM;i++){
                d=VertexCoords(j)[i]-VertexCoords(ThisVertex)[i];
                if(d>PARAMS(PolytopeEps) || d<-PARAMS(PolytopeEps)) break;
              }
              if(i==DIM){
                  report(R_err,"A vertex has been added earlier as %d\n",1+j-DIM);
                  NextVertex--;
                  return;
              }
            }
            report(R_warn,"Vertex %d is inside the approximation\n",ThisVertex-DIM+1);
        }
        return;
    }
    // some statistics
    d=(double)dd_stats.facet_pos*(double)dd_stats.facet_neg;
    if(dd_stats.max_tests<d)dd_stats.max_tests=d;
    dd_stats.avg_tests = ((dd_stats.iterations-1)*dd_stats.avg_tests+d)
          /((double)dd_stats.iterations);
    // this is the critical part: for each positive/negative facet pair
    // compute whether this is a new facet.
    // go through all negative facets
    if(DIM<=2){ // search_ridges_with() assumes DIM>=3 
        /* search_ridges() require DIM >=3. For DIM=2 f1-f2 is a ridge
           iff threre is a vertex adjacet to both of them */
        NegIdx = FacetPosnegList+(MaxFacets-1);
        for(j=0;j<dd_stats.facet_neg;j++,NegIdx--){
            if(facet_intersection(*NegIdx,*PosIdx))
                create_new_facet(*NegIdx,*PosIdx,0);
        }
    } else { search_ridges(); }
    if(OUT_OF_MEMORY || dobreak){
         dd_stats.facet_new=0;
         return;
    }
    // calculate the number of new facets to dd_stats_facet_new
    dd_stats.facet_new=NewFacet[0];
    for(i=1;i<ThreadNo;i++) dd_stats.facet_new += NewFacet[i];
    // more statistics
    i = dd_stats.facet_zero + dd_stats.facet_pos + dd_stats.facet_new;
    if(dd_stats.max_facets < i) dd_stats.max_facets = i;
    if(dd_stats.max_facetsadded<dd_stats.facet_new){
        dd_stats.max_facetsadded=dd_stats.facet_new; 
    }
    dd_stats.avg_facetsadded = ((dd_stats.iterations-1)*dd_stats.avg_facetsadded+
         dd_stats.facet_new)/((double)dd_stats.iterations);
    // delete negative facets from FacetLiving
    NegIdx = FacetPosnegList+(MaxFacets-1);
    for(j=0;j<dd_stats.facet_neg;j++,NegIdx--)
        clear_bit(FacetLiving,*NegIdx);
    // move new facets to NextFacet until there is a space
    AllNewFacets=dd_stats.facet_new;
    for(threadId=0; threadId<ThreadNo; threadId++){
        while(NewFacet[threadId]>0 && NextFacet<MaxFacets){
            NewFacet[threadId]--; AllNewFacets--;
            move_NewFacet_th(threadId,NextFacet);
            make_facet_living(NextFacet); NextFacet++;
        }
        if(NextFacet>=MaxFacets) break;
        // threadId remains where we needed it later
    }
    if(AllNewFacets==0) return; // done
    while(NewFacet[threadId]==0) threadId++;
    // no more direct space, fill holes first
    dd_stats.facets_compressed_no ++;
    // clear adjacency list of vertices
    for(i=0;i<NextVertex;i++) intersect_VertexAdj_FacetLiving(i);
    // move new facets into free facet slots
    fno=0; for(i=0;i<FacetBitmapBlockSize;i++){
        j=fno; fc=~FacetLiving[i]; // complement ...
        while(fc){
            while((fc&7)==0){ j+=3; fc>>=3; }
            if((fc&1)){
                if(j>=MaxFacets) goto finish_compress;
                NewFacet[threadId]--;
                move_NewFacet_th(threadId,j);
                make_facet_living(j);
                AllNewFacets--;
                if(AllNewFacets==0) goto finish_compress;
                while(NewFacet[threadId]==0) threadId++;
            }
            j++; fc>>=1;
        }
        fno += (1<<packshift);
    }
  finish_compress:
    if(AllNewFacets > 0){ // we still have vertices to be included
        allocate_facet_block(AllNewFacets);
        // if no memory, throw away the rest
        if(OUT_OF_MEMORY) AllNewFacets = 0;
        else while(1){
            NewFacet[threadId]--;
            move_NewFacet_th(threadId, NextFacet);
            make_facet_living(NextFacet);
            NextFacet++;
            AllNewFacets--;
            if(AllNewFacets==0) break;
            while(NewFacet[threadId]==0) threadId++;
        }
        return;
    }
    // until fno all facet slots are occupied, compress the rest
    compress_from(fno);
}


/***********************************************************************
* Interrogate facets and vertices
*
* int get_next_facet(from,to[0:DIM])
*    return the next living but not final facet number starting at
*    from. If no such vertices are, return -1.
*    When from=-1 and RandomFacet is set, pick the starting number
*    randomly. 
* bool store_vertex(double *coords)
*    called only after interrupt. Check if the vertex is new, and if yes, 
*    store it. Do not handle adjacency lists. Return 1 if this vertex
*    was not seen before; and zero otherwise.
* void make_checkpoint(void)
*    create the next checkpoint file. Fail silently
* void make_dump(void)
*    when requested by a signal, dump facets and vertices 
* int mrandom(v) -- internal
*    return a random integer between 0 and v-1 approximately uniformly. */

static inline int mrandom(int v)
{ return v<=1?0 : (random()&0x3fffffff)%v; }

int get_next_facet(int from,double *to/*[0 : DIM] */)
{int fno,i,j; BITMAP_t v;
    if(from<0){
        if(PARAMS(RandomFacet)){ // start from a random place
            fno=get_next_facet(mrandom(NextFacet),to);
            if(fno>=0) return fno;
        }
        from=0; 
    } else if(from >= NextFacet) return -1;
    dd_stats.facetenquires++;
    j=from & packmask;
    fno = from-j;
    i = fno>>packshift;
    if(j){
        v=(FacetLiving[i] & ~FacetFinal[i])>>j;
        while(v){
            while((v&7)==0){j+=3; v>>=3; }
            if(v&1){
               fno+=j;
               get_facet_into(fno,to);
               return fno;
            }
            j++; v>>=1;
        }
        fno += (1<<packshift); i++;
    }
    for(;i<FacetBitmapBlockSize;i++){
        j=0; v= FacetLiving[i] & ~FacetFinal[i];
        while(v){
            while((v&7)==0){j+=3; v>>=3; }
            if(v&1){
               fno+=j;
               get_facet_into(fno,to);
               return fno;
            }
            j++; v>>=1;
        }
        fno += (1<<packshift);
    }
    return -1; 
}

/* bool store_vertex(double *coords) */
int store_vertex(double *coords) /* store vertex if new */
{int i,j; int thisvertex; double d;
    for(i=DIM;i<NextVertex;i++){ /* check if this is a new vertex */
        for(j=0;j<DIM;j++){
            d=VertexCoords(i)[j]-coords[j];
            if(d>PARAMS(PolytopeEps) || d < -PARAMS(PolytopeEps)) break;
        }
        if(j==DIM) return 0; /* this one has been stored */
    }
    if(OUT_OF_MEMORY) return 1;
    dd_stats.vertexno++;
    if(NextVertex>=MaxVertices){
        allocate_vertex_block(0);   // do not allocate bitmaps
        if(OUT_OF_MEMORY) return 1; // out of memory
    }
    thisvertex=NextVertex; NextVertex++;
    for(i=0;i<DIM;i++){ VertexCoords(thisvertex)[i]=coords[i]; }
    return 1; /* new vertex */
}


/* void make_checkpoint(void) */
void make_checkpoint(void)
{static int version=0; // version number, increase at each iteration
    open_checkpoint(version);
    report(R_chk,"C checkpoint file #%03d for %s\n", version, PARAMS(ProblemName));
    version++;
    report(R_chk,"C vertices facets rows colums objects\n");
    report(R_chk,"N %d %d %d %d %d\n",NextVertex,get_facetnum()+1,
           PARAMS(ProblemRows),PARAMS(ProblemColumns),PARAMS(ProblemObjects));
    print_vertices(R_chk); print_facets(R_chk);
    close_checkpoint();
}

/* void make_dump(void) */
void make_dump(void)
{   open_dumpfile();
    report(R_chk,"C snapshot data for %s\n", PARAMS(ProblemName));
    report(R_chk,"C vertices faces rows colums objects\n");
    report(R_chk,"N %d %d %d %d %d\n",NextVertex,get_facetnum()+1,
           PARAMS(ProblemRows),PARAMS(ProblemColumns),PARAMS(ProblemObjects));
    print_vertices(R_chk); print_facets(R_chk);
    close_dumpfile();
}



/***********************************************************************
* Consistecy check
*
* int consistency_check() vertex bitmaps are set OK **/

int check_bitmap_consistency(void)
/* check that each facet and vertex has at least DIM adjacent items */
{int i,vno,fno,errno; int nn;
    errno=0;
    for(vno=0;vno<NextVertex;vno++){
        nn=0;
        for(i=0;i<FacetBitmapBlockSize;i++) nn+=get_bitcount(VertexAdj(vno)[i]);
        if(nn<DIM){
            report(R_err,"Vertex %d is on %d facets only\n",vno+DIM-1,nn);
            errno++;
        }
    }
    for(fno=0;fno<NextFacet;fno++) if(is_livingFacet(fno)){
        nn=0;
        for(i=0;i<VertexBitmapBlockSize;i++) nn+=get_bitcount(FacetAdj(fno)[i]);
        if(nn<DIM){
            report(R_err,"Facet %d contains %d vertices only\n",fno,nn);
            errno++;
        }
    }
    return errno;
}

static int check_facet_consistency(int fno)
// fno>=1 and living facet
{int i,vno,errno; double w, *f,*v;
    // check that all coords are >=0 and they add up to 1.0
    errno=0;
    w=-1.0; f=FacetCoords(fno);
    for(i=0;i<DIM;i++,f++){
       w += (*f);
       if(*f < PARAMS(PolytopeEps)){
           report(R_err,"Facet %d has negative coeff (%lg)\n",fno,*f);
           errno++;
       }
    }
    if(w<PARAMS(PolytopeEps) || w<-PARAMS(PolytopeEps)){
        report(R_err,"Sum of coeffs of facet %d differs from 1.0 (%lg)\n",fno,w);
        errno++;
    }
    for(vno=DIM;vno<NextVertex;vno++){
        f=FacetCoords(fno); v=VertexCoords(vno); w=0.0;
        for(i=0;i<DIM;i++,f++,v++) w += (*f)*(*v);
        w += (*f); if(w<PARAMS(PolytopeEps)){
            report(R_err,"Vertex %d is on the negative side of facet %d (%lg)\n",
                vno-DIM+1,fno,w);
            errno++;
        } else if(w<-PARAMS(PolytopeEps)){ // not adjacent
            if(extract_bit(FacetAdj(fno),vno)!=0 || 
               extract_bit(VertexAdj(vno),fno)!=0){
                  report(R_err,"Vertex %d and facet %d ar NOT adjacent (bitmap error)\n",
                     vno-DIM+1,fno);
                  errno++;
            }
        } else { // adjacent
            if(extract_bit(FacetAdj(fno),vno)==0 ||
               extract_bit(VertexAdj(vno),fno)==0){
                  report(R_err,"Vertex %d and facet %d ARE adjacent (bitmap error)\n",
                     vno-DIM+1,fno);
                  errno++;
            }
        }
    }
    return errno;
}

static void thread_check_consistency(int threadId) // o .. ThreadNo-1
{int fno,step;
    step=ThreadNo; ErrorNo[threadId]=0;
    for(fno=1+threadId;fno<NextFacet;fno+=step) if(is_livingFacet(fno))
        ErrorNo[threadId]+=check_facet_consistency(fno);
}

int check_consistency(void)
{int i,errno;
   #ifdef USETHREADS
     thread_execute(thread_check_consistency);
   #else /* ! USETHREADS */
     thread_check_consistency(0);
   #endif
     errno=check_bitmap_consistency();
     for(i=0;i<ThreadNo;i++) errno+=ErrorNo[i];
     return errno;
}

/* EOF */


