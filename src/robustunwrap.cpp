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

#include <cstdbool>
#include <cstring>

#include "raiseerror.h"
#include "multiqueue.h"
#include "robustunwrap.h"

#define PI 3.14159265358979

/*
 This code has been refactored to use a signed type for indexing with index_t.
 This is a deliberate choice following suggestion by many C++ veterans and
 mainly aims to prevent common mistakes arising from uncareful use of unsigned
 types specifcally (but not explusive) in combination with signed offsets.
 Choosing ptrdiff_t as a platform dependent signed counterpart to size_t seems
 natural, but comes with caveats:
 1. ptrdiff_t is commonly defined with the same bit-width as size_t, but this
    is not guaranteed. In case of uncommon plattforms this should be carefully
    checked to prevent overflow behaviour (if smaller than expected).
 2. Assuming same bit-width as size_t, ptrdiff_t might be only able to address
    half of the available memory (in bytes). Since it is used for indexing
    double arrays here this should not be any issue.
*/

index_t dim[3];               // size of data in x/y/z direction
index_t sze;                  // size of data as flattened array
const double *phase, *mag;    // input data phase/magnitude
double *unwrapped;            // global results: unwrapped phase (memory allocated by caller)
bool *flag;                   // indicator of points processed
index_t m_bsx, m_bsy, m_bsz;  // stride for x/y/z dimension in flattened data array

struct QUEUEENTRY {  // describes one voxel of the data
    index_t x, y, z; // index in 3D array
    index_t p;       // index in flattened array
    double v;        // value
};

/*
 perform unwrapping step on neighbor of qe at offset (offp, offx, offy, offz)
 checks are performed if offset points to a valid position and if data at offset
 has already been processed (flag(p) == true). In both cases the function returns
 without any action. If checks pass a potential shift by multiples of 2*PI is
 estimated and applied based on the value difference of qe and the neighbor (the
 multiple might be 0). The resulting entry is written to unwrapped (the globally
 defind output of the algorithm) and pushed to queue primaryqueuenum.

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
void Check(queue_index_t primaryqueuenum, QUEUEENTRY *qe, index_t offp, index_t offx, index_t offy, index_t offz) {
    /* first check bounds */
    QUEUEENTRY nqe{};
    nqe.x = qe->x + offx;
    if ((nqe.x < 0) || (nqe.x >= dim[0])) return;
    nqe.y = qe->y + offy;
    if ((nqe.y < 0) || (nqe.y >= dim[1])) return;
    nqe.z = qe->z + offz;
    if ((nqe.z < 0) || (nqe.z >= dim[2])) return;

    nqe.p = qe->p + offp;

    if (flag[nqe.p]) return; /* Already been here */

    /* Actually do unwrap */
    int wholepis;
    wholepis = int((phase[nqe.p] - qe->v) / PI);

    if (wholepis >= 1)
        nqe.v = (double) phase[nqe.p] - 2 * PI * int((wholepis + 1) / 2);
    else if (wholepis <= -1)
        nqe.v = (double) phase[nqe.p] + 2 * PI * int((1 - wholepis) / 2);
    else nqe.v = phase[nqe.p];

    unwrapped[nqe.p] = nqe.v;
    flag[nqe.p] = true;

    if (Push(primaryqueuenum, &nqe))
        raiseerror("Out of memory!");
}

/*
 The main unwrapping algorithm
 1. Define the polefield thresholds based on inverse magnitude image and num_unwrapbins
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
 num_unwrapbins: number of bins to use for this algorithm. this determines the
                 granularity of the polefield thresholds (must be positive)

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
void unwrap(index_t seedx, index_t seedy, index_t seedz, index_t num_unwrapbins) {
    queue_index_t unwrapbins = static_cast<queue_index_t>(num_unwrapbins);
    if (unwrapbins > NUMQUEUES) {
        raiseerror("The maximum number of unwrapbins is hardcoded to a smaller number.");
        return;
    }
    /* Minimum number of unwrapping bins is 2 */
    if (unwrapbins < 2) unwrapbins = 2;
    /* Find min and max */
    double min, max;
    min = (double) 1e38;
    max = 0;

    for (index_t i = 0; i < sze; i++) {
        min = mag[i] < min ? mag[i] : min;
        max = mag[i] > max ? mag[i] : max;
    }

    double diff = (double) 1.00001 * (max - min);

    index_t seedp = seedx + dim[0] * (seedy + dim[1] * seedz);

    flag = new bool[sze];

    if (!flag) raiseerror("Out of memory. ");
    for (index_t i = 0; i < sze; i++)
        flag[i] = false;

    for (queue_index_t i = 0; i < unwrapbins; i++)
        InitQueue(i, sizeof(QUEUEENTRY));
    QUEUEENTRY qe{};

    qe.p = seedp;
    qe.x = seedx;
    qe.y = seedy;
    qe.z = seedz;
    unwrapped[seedp] = qe.v = phase[seedp];
    flag[seedp] = true;

/* push seed */
    Push(0, &qe);

/* First, work out pole field threshold that we're going to use */
    auto *polefieldthresholds = new double[static_cast<size_t>(unwrapbins)];

    for (queue_index_t i = 0; i < unwrapbins; i++)
        polefieldthresholds[i] = min + (diff * static_cast<double>(i) / (static_cast<double>(unwrapbins) - 1.0));

    // TODO: Display warning if seedpoint is chosen badly (i.e. mag[seed]>polefieldthresholds[0])

    for (queue_index_t i = 0; i < unwrapbins; i++) {
        while (!IsEmpty(i)) {
            Pop(i, &qe);
            if (mag[qe.p] >
                polefieldthresholds[i]) { /* too close to a scary pole, so just defer by pushing to other stack */
                queue_index_t ind;
                for (ind = i + 1; mag[qe.p] > polefieldthresholds[ind]; ind++);
                Push(ind, &qe);                                        /* just defer by pushing to relevant stack */
            } else {
                Check(i, &qe, +m_bsz, 0, 0, 1);
                Check(i, &qe, -m_bsz, 0, 0, -1);
                Check(i, &qe, +m_bsy, 0, 1, 0);
                Check(i, &qe, -m_bsy, 0, -1, 0);
                Check(i, &qe, +m_bsx, 1, 0, 0);
                Check(i, &qe, -m_bsx, -1, 0, 0);
            }
        }
        TerminateQueue(i);    /* done with this Queue so free up memory */
    }

    delete flag;
}

/*
 Wrapper to initialize the global variables used in this source file while
 exposing a complete function API to call the unwrapping algorithm.

 TODO: This function should not be needed and disappear after thorough refactoring.
*/
void unwrap_helper(index_t seedx,
                   index_t seedy,
                   index_t seedz,
                   index_t num_unwrapbins,
                   const index_t *dim_,
                   index_t sze_,
                   index_t m_bsx_,
                   index_t m_bsy_,
                   index_t m_bsz_,
                   const double *phase_,
                   const double *mag_,
                   double *unwrapped_) {
    sze = sze_;
    m_bsx = m_bsx_;
    m_bsy = m_bsy_;
    m_bsz = m_bsz_;
    phase = phase_;
    mag = mag_;
    unwrapped = unwrapped_; // memory allocated/managed by calling interface
    {
        size_t dim_sz = sizeof(dim) / sizeof(dim[0]);
        for (size_t i = 0; i < dim_sz; i++)
            dim[i] = dim_[i];
    }
    //memcpy(x,m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],m_chunk[queuenum]);

    unwrap(seedx, seedy, seedz, num_unwrapbins);
}
