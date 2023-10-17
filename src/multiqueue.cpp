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
#include <new>
#include "raiseerror.h"
#include "multiqueue.h"

static queue_index_t m_bot[NUMQUEUES], m_top[NUMQUEUES]; // FI and LI element for each queue
static queue_index_t m_size[NUMQUEUES]; // memory size of each queue in elements
static size_t m_sizeb[NUMQUEUES], m_chunk[NUMQUEUES]; // memory size of each queue in bytes (m_sizeb) and size of a single element (m_chunk) for each queue in bytes
static char *m_q[NUMQUEUES]; // pointer to allocated memory for each queue (or 0)

/*
 Initialize the queue at index queuenum for elements of size chunk bytes
 This function only sets the metadata, the actual memory allocation is done
 in lazy manner by push()
*/
void InitQueue(queue_index_t queuenum, size_t chunk) {
    m_chunk[queuenum] = chunk;
    m_bot[queuenum] = 0;
    m_top[queuenum] = 0;
    m_size[queuenum] = DEFAULTBLOCK;
    m_sizeb[queuenum] = m_size[queuenum] * m_chunk[queuenum];
    m_q[queuenum] = nullptr;
}

/*
 free the memory for queue at index queuenum and re-initialize as emty queue
*/
void TerminateQueue(queue_index_t queuenum) {
    if (m_q[queuenum]) {
        delete m_q[queuenum];
        InitQueue(queuenum, m_chunk[queuenum]);
    }
}

/*
 Initialize the first num_queues for elements of size chunk bytes
 This is a wrapper around InitQueue for conveninece

 returns 1 on success
 returns 0 if num_queues is larger than NUMQUEUES
*/
int InitNQueues(queue_index_t num_queues, size_t chunk) {
    if (num_queues > NUMQUEUES)
        return 0;
    for (queue_index_t i = 0; i < num_queues; i++)
        InitQueue(i, chunk);
    return 1;
}

/*
 push element x into queue at index queuenum
 this function takes care of allocating memory and expanding the queue size
 if necessary.
 x is copied from the given pointer. size is assumed to be m_chunk[queuenum]
 bytes as defined with InitQueue()
*/
int Push(queue_index_t queuenum, void *x) {
    if (!m_q[queuenum]) m_q[queuenum] = new(std::nothrow) char[m_sizeb[queuenum]];

    // TODO: Consider getting rid of nothrow and pass the exception
    if (!m_q[queuenum]) {
        raiseerror("Out of memory - could not generate new point queue.");
        return (-1); // should never be reached
    }

    memcpy(m_q[queuenum] + m_top[queuenum] * m_chunk[queuenum], x, m_chunk[queuenum]);
    m_top[queuenum]++;

    if (m_top[queuenum] == m_size[queuenum]) m_top[queuenum] = 0;

    if (m_top[queuenum] == m_bot[queuenum]) {
        char *newq;
        size_t newsize, newsizeb, abovebot_b;
        // previously out of memory - now auto-expand
        newsize = m_size[queuenum] + BLOCKINCREMENT;
        newsizeb = newsize * m_chunk[queuenum];
        newq = new(std::nothrow) char[newsizeb];

        if (!newq) {
            raiseerror("Out of memory - point queue full.");
            return (-1);
        }

        // While we're shifting circular buffer, it is actually easier to re-origin it
        // to zero.
        // first, copy top bit from m_bot upwards
        abovebot_b = (m_size[queuenum] - m_bot[queuenum]) * m_chunk[queuenum];
        if (abovebot_b != 0) memcpy(newq, (char *) m_q[queuenum] + m_bot[queuenum] * m_chunk[queuenum], abovebot_b);
        // then, do m_top downwards
        if (m_top[queuenum] != 0)
            memcpy((char *) newq + abovebot_b, m_q[queuenum], m_top[queuenum] * m_chunk[queuenum]);
        m_bot[queuenum] = 0;
        m_top[queuenum] = m_size[queuenum];
        m_size[queuenum] = newsize;
        m_sizeb[queuenum] = newsizeb;
        delete[] m_q[queuenum]; // recover old memory
        m_q[queuenum] = newq;
    }

    return (0);
}

/*
 Copy first element of the queue to pointer x and remove from the queue
 returns 0 on success
 returns 1 if queue is empty (and raises an error)
*/
int Pop(queue_index_t queuenum, void *x) {
    if (m_bot[queuenum] == m_top[queuenum]) {
        raiseerror("Out of bounds - trying to Pop() empty queue.");
        return (1);
    }

    memcpy(x, m_q[queuenum] + m_bot[queuenum] * m_chunk[queuenum], m_chunk[queuenum]);
    m_bot[queuenum]++;
    if (m_bot[queuenum] == m_size[queuenum]) m_bot[queuenum] = 0;

    return (0);
}

/*
 Check if queue at index is empty
 This function assumes the queue is initialized.
 If InitQueue() has not been caled for queuenum behaviour is undefined.

 returns 1 for empty queue
 retruns 0 otherwise
*/
int IsEmpty(queue_index_t queuenum) {
    if (m_bot[queuenum] == m_top[queuenum])
        return (1);
    else
        return (0);
}