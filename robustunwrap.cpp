/*
 3D Unwrapping algorithm
 Rhodri Cusack 2000-2006
 Algorithm described in 
 Cusack, R. & Papadakis, N. (2002) 
    "New robust 3-D phase unwrapping algorithms: application to magnetic field mapping and undistorting echoplanar images."
    Neuroimage. 2002 Jul;16(3 Pt 1):754-64.
 Distributed under MIT License
 Comments to cusackrh@tcd.ie
 Version 3.00 Oct 2006 adapted for matlab
 */

#include <string.h> // memcpy

#include "raiseerror.h"
#include "robustunwrap.h"

#define PI 3.14159265358979
#define FILEFORMAT_CAMRES 1
#define FILEFORMAT_ANALYZE 2
#define MAXFILENAMESIZE 1024

int errcode;            // indicator of an error
ptrdiff_t dim[3];       // size of data in x/y/z direction
ptrdiff_t sze;             // size of data as fattened array
const double *phase,*mag;     // input data phase/magnitude
double *unwrapped;      // global results: unwrapped phase
bool *flag;             // indicator of points processed
ptrdiff_t m_bsx,m_bsy,m_bsz; // stride for x/y/z dimension in flattened data array

struct QUEUEENTRY {
    ptrdiff_t x,y,z;
    ptrdiff_t p;
    double v;
};

/*
 wrapper to raise an errors in API specific ways
 raiseerror_api imported from raiseerror.h should be linked to raise errors
 in a way specific to the frontend (e.g. Matlab's mexErrMsgTxt)

 Direct global variable access:
 errorcode (indicator of an error): is set
*/
void raiseerror(const char *msg) {
    raiseerror_api(msg);
    errcode=1;
};

/* QUEUE
 The follwoing functions implement an array of NUMQUEUES FIFO queues for
 QUEUEENTRY elements. The individual queues get initialized with space for
 DEFAULTBLOCK number of elements and grow by BLOCKINCREMENT if full
 automatically.
 The initial potion of memeory is reserved in lazy manner upon first push to the spcific queue
*/
const long NUMQUEUES=10000; // number of independent queues available
const long DEFAULTBLOCK=100; // initial size of each queue in elements
const long BLOCKINCREMENT=500; // automatic expansion size if queue is full

long m_bot[NUMQUEUES], m_top[NUMQUEUES]; // FI and LI element for each queue
long m_size[NUMQUEUES], m_chunk[NUMQUEUES], m_sizeb[NUMQUEUES]; // memory size of each queue in elements (m_size) and bytes (m_sizeb) and size of a single element (m_chunk) for each queue
char *m_q[NUMQUEUES]; // pointer to allocated memory for each queue (or 0)

/*
 Initialize the queue at index queuenum for elements of size chunk bytes
 This function only sets the metadata, the actual memory allocation is done
 in lazy manner by push()
*/
void InitQueue(int queuenum,int chunk) {
    m_chunk[queuenum]=chunk;
    m_bot[queuenum]=0;
    m_top[queuenum]=0;
    m_size[queuenum]=DEFAULTBLOCK;
    m_sizeb[queuenum]=m_size[queuenum]*m_chunk[queuenum];
    m_q[queuenum]=0;
};

/*
 free the memory for queue at ince queuenum
 other metadata like m_bot and m_top are not reset so using Pop() or Push()
 for queuenum after calling this function is not safe.
*/
void TerminateQueue(int queuenum) {
    if (m_q[queuenum]) delete m_q[queuenum];
    m_q[queuenum]=0;
};

/*
 push element x into queue at index queuenum
 this function takes care of allocating memory and expanding the queue size
 if necessary.
 x is copied from the given pointer. size is assumed to be m_chunk[queuenum]
 bytes as defined with InitQueue()
*/
int Push(int queuenum,void *x) {
    if (!m_q[queuenum]) m_q[queuenum]=new char[m_sizeb[queuenum]];
    if (!m_q[queuenum]) {
        raiseerror("Out of memory - could not generate new point queue.");
        return(-1); // should never be reached
    };

    int stacksize;
    stacksize=m_top[queuenum]-m_bot[queuenum];
    if (stacksize<0) stacksize+=m_size[queuenum];
    memcpy(m_q[queuenum]+m_top[queuenum]*m_chunk[queuenum],x,m_chunk[queuenum]);
    m_top[queuenum]++;

    if (m_top[queuenum]==m_size[queuenum]) m_top[queuenum]=0;

    if (m_top[queuenum]==m_bot[queuenum]) {
        char *newq;
        long newsize, newsizeb, abovebot_b;
        // previously out of memory - now auto-expand
        newsize=m_size[queuenum]+BLOCKINCREMENT;
        newsizeb=newsize*m_chunk[queuenum];
        newq=new char[newsizeb];
        if (!newq) {
            raiseerror("Out of memory - point queue full.");
            return(-1); // should be never reached
        };
        // While we're shifting circular buffer, it is actually easier to re-origin it
        // to zero.
        // first, copy top bit from m_bot upwards
        abovebot_b=(m_size[queuenum]-m_bot[queuenum])*m_chunk[queuenum];
        if (abovebot_b>0) memcpy(newq,(char *)m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],abovebot_b);
        // then, do m_top downwards
        if (m_top[queuenum]!=0) memcpy((char *) newq+abovebot_b,m_q[queuenum],m_top[queuenum]*m_chunk[queuenum]);
        m_bot[queuenum]=0;
        m_top[queuenum]=m_size[queuenum];
        m_size[queuenum]=newsize;
        m_sizeb[queuenum]=newsizeb;
        delete m_q[queuenum]; // recover old memory
        m_q[queuenum]=newq;
    };

    return(0);
};

/*
 Copy first element of the queue to pointer x and remove from the queue
 returns 0 on success
 returns 1 if queue is empty
*/
int Pop(int queuenum,void *x) {
    if (m_bot[queuenum]==m_top[queuenum]) return(1);

    memcpy(x,m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],m_chunk[queuenum]);
    m_bot[queuenum]++;
    if (m_bot[queuenum]==m_size[queuenum]) m_bot[queuenum]=0;

    return(0);
};
/* QUEUE end */

/*
 perform unwrapping step on neighbor of qe at offset (offp, offx, offy, offz)
 checks are performed if offset points to a valid positon amd if data at offset
 has already been processed (flag(p) == true). In both cases the function returns
 without any action. If checks oass a potential shift by multiples of 2*PI is
 estimated and applied based on the value difference of qe and the neighbor (the
 multiple might be 0). The resulting entry is written to unwrapped (the globally
 defind output of the algorithm) and and pushed to queue primaryqueuenum.

 Arguments:
 primaryqueuenum: index of the queue for pushing the resulting entry
 qe: entry of the current position
 offp: offset to the current position as flatened array index (must correspond to offx/y/z)
 offx/y/z: offset to the current position along x/y/z (or dim[0]/[1]/[2]) axis (must correspond to offp)

 Direct global variable access:
 unwrapped (global results): unwrapped value for point at offset is written
 phase (the input phase data): read
 dim (size of data in x/y/z direction): read
 flag (indicator of points processed): flag for point at offset is set
*/
void Check(int primaryqueuenum, QUEUEENTRY* qe,long offp, int offx, int offy, int offz) {
    /* first check bounds */
    QUEUEENTRY nqe;
    nqe.x=qe->x+offx;
    if ((nqe.x<0)||(nqe.x>=dim[0])) return;
    nqe.y=qe->y+offy;
    if ((nqe.y<0)||(nqe.y>=dim[1])) return;
    nqe.z=qe->z+offz;
    if ((nqe.z<0)||(nqe.z>=dim[2])) return;

    nqe.p=qe->p+offp;

    if (flag[nqe.p]) return; /* Already been here */

    /* Actually do unwrap */
    int wholepis;
    wholepis=int((phase[nqe.p]-qe->v)/PI);

    if (wholepis>=1)
        nqe.v=(double) phase[nqe.p]-2*PI*int((wholepis+1)/2);
    else if (wholepis<=-1)
        nqe.v=(double) phase[nqe.p]+2*PI*int((1-wholepis)/2);
    else nqe.v=phase[nqe.p];

    unwrapped[nqe.p]=nqe.v;
    flag[nqe.p]=true;

    if (Push(primaryqueuenum,&nqe))
        raiseerror("Out of memory!");
};

/*
 The main unwrapping algorithm
 1. Define the polefield thresholds based on inverse magnitude image and UNWRAPBINS
 2. Initialize flag and queues (and queue seedpoint for processing)
 3. start processing by popping an entry (current point) from the lowest indexed
    queue and
     a. put into a higher indexed queue if value is above polefield threshold or
     b. unwrap and queue the 6 closest neighbors otherwise
 4. once the queue is empty terminate it, increment the queue index and repeat 3.
 5. if all queues have been emptied the unwrapping is done.

 The given seed is not checked for polefield threshold. This might result in jumping
 into a high unwrapbin right away if seed is close to a pole and ultimately lead to
 underperformance.

 Arguments:
 seedx/y/z: coordinates (0-based indices) of the seedpoint along dim[0]/[1]/[2]
 UNWRAPBINS: number of bins to use for this algorithm. this determines the granularity
             of the polefield thresholds

 Direct global variable access:
 unwrapped (global results): seedpoint is written
 phase (the input phase data): seedpoint is read
 mag (the input magnitude data): read
 dim (size of data in x/y/z direction): read
 sze (size of data as fattened array): read
 m_bsx/m_bsy/m_bsz (stride for x/y/z dimension): read
 flag (indicator of points processed): is initialized
 errorcode (indicator of an error): initialize and read
*/
void unwrap(ptrdiff_t seedx, ptrdiff_t seedy, ptrdiff_t seedz, ptrdiff_t UNWRAPBINS) {
    ptrdiff_t i;
    ptrdiff_t seedp;

    errcode=0;

    /* Minimum number of unwrapping bins is 2 */
    if (UNWRAPBINS<2) UNWRAPBINS=2;
    /* Find min and max */
    double min,max;
    min=(double) 1e38;
    max=0;

    for (i=0; i<sze; i++) {
        min=mag[i]<min?mag[i]:min;
        max=mag[i]>max?mag[i]:max;
    };

    double diff = (double) 1.00001 * (max - min);

    seedp=seedx+dim[0]*(seedy+dim[1]*seedz);

    flag=new bool[sze];

    if (!flag) raiseerror("Out of memory. ");
    for (i=0; i<sze; i++)
        flag[i]=false;

    for (i=0; i<UNWRAPBINS; i++)
        InitQueue(i,sizeof(QUEUEENTRY));
    QUEUEENTRY qe;

    qe.p=seedp;
    qe.x=seedx;
    qe.y=seedy;
    qe.z=seedz;
    unwrapped[seedp]=qe.v=phase[seedp];
    flag[seedp]=true;

/* push seed */
    Push(0,&qe);

/* First, work out pole field threshold that we're going to use */
    double *polefieldthresholds = new double [UNWRAPBINS];

    for (i=0; i<(UNWRAPBINS); i++)
        polefieldthresholds[i] = min + (diff * i / (UNWRAPBINS - 1)); 

    // TODO: Display warning if seedpoint is chosen badly (i.e. mag[seed]>polefieldthresholds[0])

    for (i=0; i<UNWRAPBINS ; i++) {
        while (!errcode && !Pop(i,&qe)) {
            if (mag[qe.p]>polefieldthresholds[i]) { /* too close to a scary pole, so just defer by pushing to other stack */
                int ind;
                for (ind = i+1; mag[qe.p] > polefieldthresholds[ind]; ind++);
                Push(ind, &qe);										/* just defer by pushing to relevant stack */
            } else {
                Check(i,&qe,+m_bsz,0,0,1);
                Check(i,&qe,-m_bsz,0,0,-1);
                Check(i,&qe,+m_bsy,0,1,0);
                Check(i,&qe,-m_bsy,0,-1,0);
                Check(i,&qe,+m_bsx,1,0,0);
                Check(i,&qe,-m_bsx,-1,0,0);
            };
        };
        TerminateQueue(i);	/* done with this Queue so free up memory */
        
        if (errcode) break;
    };
    
    delete flag;
};

/*
 Wrapper to initialize the global variables used in this source file while
 exposing a complete function API to call the unwrapping algorithm.

 TODO: This function should not be needed and disappear after thorough refactoring.
*/
void unwrap_helper(const ptrdiff_t seedx, const ptrdiff_t seedy, const ptrdiff_t seedz,
                   const ptrdiff_t UNWRAPBINS,
                   const ptrdiff_t *dim_,
                   const ptrdiff_t sze_,
                   const ptrdiff_t m_bsx_,
                   const ptrdiff_t m_bsy_,
                   const ptrdiff_t m_bsz_,
                   const double *phase_,
                   const double *mag_,
                   double *unwrapped_) {
    sze = sze_;
    m_bsx = m_bsx_;
    m_bsy = m_bsy_;
    m_bsz = m_bsz_;
    phase = phase_;
    mag = mag_;
    unwrapped = unwrapped_;
    {
        size_t dim_sz = sizeof(dim) / sizeof(dim[0]);
        for (size_t i = 0; i < dim_sz; i++)
            dim[i] = dim_[i];
    }
    //memcpy(x,m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],m_chunk[queuenum]);

    if (UNWRAPBINS > NUMQUEUES) {
        raiseerror("The maximum number of UNWRAPBINS is hardcoded to a smaller number.");
        return;
    }

    unwrap(seedx, seedy, seedz, UNWRAPBINS);
};
