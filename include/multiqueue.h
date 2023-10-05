/* MULTIQUEUE
 implement an array of NUMQUEUES FIFO queues for
 QUEUEENTRY elements. The individual queues get initialized with space for
 DEFAULTBLOCK number of elements and grow by BLOCKINCREMENT if full
 automatically.
 The initial potion of memory is reserved in lazy manner upon first push to the specific queue
*/

#pragma once

const long NUMQUEUES = 10000; // number of independent queues available
const long DEFAULTBLOCK = 100; // initial size of each queue in elements
const long BLOCKINCREMENT = 500; // automatic expansion size if queue is full

/*
 Initialize the queue at index queuenum for elements of size chunk bytes
 This function only sets the metadata, the actual memory allocation is done
 in lazy manner by push()
*/
void InitQueue(int queuenum, int chunk);

/*
 free the memory for queue at index queuenum
 other metadata like m_bot and m_top are not reset so using Pop() or Push()
 for queuenum after calling this function is not safe.
*/
void TerminateQueue(int queuenum);

/*
 push element x into queue at index queuenum
 this function takes care of allocating memory and expanding the queue size
 if necessary.
 x is copied from the given pointer. size is assumed to be m_chunk[queuenum]
 bytes as defined with InitQueue()
*/
int Push(int queuenum, void *x);

/*
 Copy first element of the queue to pointer x and remove from the queue
 returns 0 on success
 returns 1 if queue is empty
*/
int Pop(int queuenum, void *x);
