/* MULTIQUEUE
 implement an array of NUMQUEUES FIFO queues for
 QUEUEENTRY elements. The individual queues get initialized with space for
 DEFAULTBLOCK number of elements and grow by BLOCKINCREMENT if full
 automatically.
 The initial potion of memory is reserved in lazy manner upon first push to the specific queue
*/

#pragma once

typedef size_t queue_index_t; // type used for indexing and counting.
                              // this is assumed to be compatible with size_t

const queue_index_t NUMQUEUES = 10000; // number of independent queues available
const queue_index_t DEFAULTBLOCK = 100; // initial size of each queue in elements
const queue_index_t BLOCKINCREMENT = 500; // automatic expansion size in elements if queue is full

/*
 Initialize the queue at index queuenum for elements of size chunk bytes
 This function only sets the metadata, the actual memory allocation is done
 in lazy manner by push()
*/
void InitQueue(queue_index_t queuenum, size_t chunk);

/*
 free the memory for queue at index queuenum
 other metadata like m_bot and m_top are not reset so using Pop() or Push()
 for queuenum after calling this function is not safe.
*/
void TerminateQueue(queue_index_t queuenum);

/*
 push element x into queue at index queuenum
 this function takes care of allocating memory and expanding the queue size
 if necessary.
 x is copied from the given pointer. size is assumed to be m_chunk[queuenum]
 bytes as defined with InitQueue()
*/
int Push(queue_index_t queuenum, void *x);

/*
 Copy first element of the queue to pointer x and remove from the queue
 returns 0 on success
 returns 1 if queue is empty (and raises an error)
*/
int Pop(queue_index_t queuenum, void *x);

/*
 Check if queue at index is empty
 This function assumes the queue is initialized.
 If InitQueue() has not been caled for queuenum behaviour is undefined.

 returns 1 for empty queue
 retruns 0 otherwise
*/
int IsEmpty(queue_index_t queuenum);
