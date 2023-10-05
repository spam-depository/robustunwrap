/* MULTIQUEUE
 The following functions implement an array of NUMQUEUES FIFO queues for
 QUEUEENTRY elements. The individual queues get initialized with space for
 DEFAULTBLOCK number of elements and grow by BLOCKINCREMENT if full
 automatically.
 The initial potion of memory is reserved in lazy manner upon first push to the specific queue

 Extracted from robustunwrap
 (c) Rhodri Cusack 2000-2006
 Comments to cusackrh@tcd.ie
 Distributed under MIT License
*/

#include <cstddef>
#include <cstring>
#include "raiseerror.h"
#include "multiqueue.h"

static long m_bot[NUMQUEUES], m_top[NUMQUEUES]; // FI and LI element for each queue
static long m_size[NUMQUEUES], m_chunk[NUMQUEUES], m_sizeb[NUMQUEUES]; // memory size of each queue in elements (m_size) and bytes (m_sizeb) and size of a single element (m_chunk) for each queue
static char *m_q[NUMQUEUES]; // pointer to allocated memory for each queue (or 0)

/*
 Initialize the queue at index queuenum for elements of size chunk bytes
 This function only sets the metadata, the actual memory allocation is done
 in lazy manner by push()
*/
void InitQueue(int queuenum, int chunk) {
    m_chunk[queuenum] = chunk;
    m_bot[queuenum] = 0;
    m_top[queuenum] = 0;
    m_size[queuenum] = DEFAULTBLOCK;
    m_sizeb[queuenum] = m_size[queuenum] * m_chunk[queuenum];
    m_q[queuenum] = nullptr;
}

/*
 free the memory for queue at ince queuenum
 other metadata like m_bot and m_top are not reset so using Pop() or Push()
 for queuenum after calling this function is not safe.
*/
void TerminateQueue(int queuenum) {
    if (m_q[queuenum]) delete m_q[queuenum];
    m_q[queuenum] = nullptr;
}

/*
 push element x into queue at index queuenum
 this function takes care of allocating memory and expanding the queue size
 if necessary.
 x is copied from the given pointer. size is assumed to be m_chunk[queuenum]
 bytes as defined with InitQueue()
*/
int Push(int queuenum, void *x) {
    if (!m_q[queuenum]) m_q[queuenum] = new char[m_sizeb[queuenum]];

    // TODO: Can be removed. "new char[..]" will throw an exception when out of memory. It will not return null_ptr.
    // See https://en.cppreference.com/w/cpp/memory/new/operator_new
    if (!m_q[queuenum]) {
        raiseerror("Out of memory - could not generate new point queue.");
        return (-1); // should never be reached
    }

    // TODO: What exactly is stacksize supposed to do? It's never used anywhere else although it is assigned here.
    long stacksize;
    stacksize = m_top[queuenum] - m_bot[queuenum];
    if (stacksize < 0) stacksize += m_size[queuenum];

    memcpy(m_q[queuenum] + m_top[queuenum] * m_chunk[queuenum], x, m_chunk[queuenum]);
    m_top[queuenum]++;

    if (m_top[queuenum] == m_size[queuenum]) m_top[queuenum] = 0;

    if (m_top[queuenum] == m_bot[queuenum]) {
        char *newq;
        long newsize, newsizeb, abovebot_b;
        // previously out of memory - now auto-expand
        newsize = m_size[queuenum] + BLOCKINCREMENT;
        newsizeb = newsize * m_chunk[queuenum];
        newq = new char[newsizeb];

        // TODO: Can be removed. "new char[..]" will throw an exception when out of memory. It will not return null_ptr.
        if (!newq) {
            raiseerror("Out of memory - point queue full.");
            return (-1); // should be never reached
        }

        // While we're shifting circular buffer, it is actually easier to re-origin it
        // to zero.
        // first, copy top bit from m_bot upwards
        abovebot_b = (m_size[queuenum] - m_bot[queuenum]) * m_chunk[queuenum];
        if (abovebot_b > 0) memcpy(newq, (char *) m_q[queuenum] + m_bot[queuenum] * m_chunk[queuenum], abovebot_b);
        // then, do m_top downwards
        if (m_top[queuenum] != 0)
            memcpy((char *) newq + abovebot_b, m_q[queuenum], m_top[queuenum] * m_chunk[queuenum]);
        m_bot[queuenum] = 0;
        m_top[queuenum] = m_size[queuenum];
        m_size[queuenum] = newsize;
        m_sizeb[queuenum] = newsizeb;
        // TODO: Seems wrong as m_q[index] is an array. Use delete[]
        delete[] m_q[queuenum]; // recover old memory
        m_q[queuenum] = newq;
    }

    return (0);
}

/*
 Copy first element of the queue to pointer x and remove from the queue
 returns 0 on success
 returns 1 if queue is empty
*/
int Pop(int queuenum, void *x) {
    if (m_bot[queuenum] == m_top[queuenum]) return (1);

    memcpy(x, m_q[queuenum] + m_bot[queuenum] * m_chunk[queuenum], m_chunk[queuenum]);
    m_bot[queuenum]++;
    if (m_bot[queuenum] == m_size[queuenum]) m_bot[queuenum] = 0;

    return (0);
}
