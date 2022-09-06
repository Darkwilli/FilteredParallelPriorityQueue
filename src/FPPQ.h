#pragma once

#include <atomic>
#include <vector>
#include <execution>
#include <algorithm>

#include <cassert>

#include <omp.h>
#include <Windows.Foundation.h>
/*
* 
* A General Implementation of the Filtered Parallel Priority Queue 
* 
* 
*/

template<typename T>
class FPPQ {
protected:
	const int maxPopIndex = 31; // Is actually the count of nodes that are allowed to be popped
	const int filterLevelCount = 4;
public:
	FPPQ(int maxCapacity);
	// Fills The PQ with a presorted array does use more threads than countThreads because of the parallel std::sort
	// To make accurate benchmarks use a custom parallel sorting implementation instead
	void fill(std::vector<std::pair<T, int>> errorValues, int countThreads);

	virtual int pop(); // pops a top element

	// Updates the element with the id: "lookupId" with a new errorValue. This id can be used to get the current node of the element, if it is currently updated in the upward direction (see original paper) and the current error.
	// fromInsertion marks if the update comes from an insert operation, which allows to skip a few steps
	virtual void update(int lookupId, T errorValue, bool fromInsertion = false);

	// Insert the element with the id: "lookupId" with an error of errorValue.
	virtual void insert(int lookupId, T errorValue);

	// Locks the node with id: "nodeId"
	inline void lock(int nodeId) {
		EnterCriticalSection(&m_locks[correctLockId(nodeId)]);
		//while(!tryLock(nodeId))
			//m_locks[correctNodeId(nodeId)].wait(true);
		
	}

	// Unlocks the node with id: "nodeId"
	inline void unlock(int nodeId) {
		LeaveCriticalSection(&m_locks[correctLockId(nodeId)]);
		//m_locks[correctNodeId(nodeId)].clear(std::memory_order_release);
		//m_locks[correctNodeId(nodeId)].notify_one();
	}

	// Trys to lock the node with id: "nodeId" and returns if successfull
	inline bool tryLock(int nodeId) {
		return TryEnterCriticalSection(&m_locks[correctLockId(nodeId)]);
		//return !m_locks[correctNodeId(nodeId)].test_and_set(std::memory_order_acquire);
	}

	inline int getSize()
	{
		return static_cast<int>(m_nodes.size());
	}
	inline std::vector<int>& getNodes()
	{
		return m_nodes;
	}

	// Checks if the indexing is still valid (if the child of a parent has the parent as a parent etc.)
	void debugCheckIndices() {
		for(int i = 0; i < filteredElements - filterLevelSize; i++) {
			int childNode1 = leftChildNodeIdAboveFilter(i);
			int childNode2 = rightChildNodeIdAboveFilter(i);
			int parentNode1 = parentNodeIdAboveFilter(childNode1);
			int parentNode2 = parentNodeIdAboveFilter(childNode2);
			if(parentNode1 != parentNode2 || parentNode1 != i)
				__debugbreak();
		}
		for(int i = filteredElements - filterLevelSize; i < filterLevelEnd - filterLevelSize; i++) {
			volatile int childNode1, childNode2, filterLevelOffset;
			filterLevelOffset = filterLevelDown(i);
			if(filterConditionDown(i, filterLevelOffset))
			{
				childNode2 = leftChildNodeId2InFilter(i, filterLevelOffset);
				childNode1 = childNodeId1InFilter(i);
			} else {
				childNode1 = childNodeId1InFilter(i);
				childNode2 = rightChildNodeId2InFilter(i, filterLevelOffset);
			}
			{
				volatile int leftParentNodeId, rightParentNodeId;
				int nodeId = childNode1;
				filterLevelOffset = filterLevelUp(nodeId);
				if(filterConditionUp(nodeId, filterLevelOffset))
				{
					leftParentNodeId = parentNodeId1InFilter(nodeId);
					rightParentNodeId = rightParentNodeId2InFilter(nodeId, filterLevelOffset);
				} else {
					leftParentNodeId = leftParentNodeId2InFilter(nodeId, filterLevelOffset);
					rightParentNodeId = parentNodeId1InFilter(nodeId);
				}
				if(!(leftParentNodeId == i || rightParentNodeId == i))
					__debugbreak();
			}
			{
				volatile int leftParentNodeId, rightParentNodeId;
				int nodeId = childNode2;
				filterLevelOffset = filterLevelUp(nodeId);
				if(filterConditionUp(nodeId, filterLevelOffset))
				{
					leftParentNodeId = parentNodeId1InFilter(nodeId);
					rightParentNodeId = rightParentNodeId2InFilter(nodeId, filterLevelOffset);
				} else {
					leftParentNodeId = leftParentNodeId2InFilter(nodeId, filterLevelOffset);
					rightParentNodeId = parentNodeId1InFilter(nodeId);
				}
				if(!(leftParentNodeId == i || rightParentNodeId == i))
					__debugbreak();
			}


		}
		for(int i = filterLevelEnd - filterLevelSize; i < 1000000; i++) {
			if(i < 0)
				i = 0;
			int childNode1 = leftChildNodeIdBelowFilter(i);
			int childNode2 = rightChildNodeIdBelowFilter(i);
			int parentNode1 = parentNodeIdBelowFilter(childNode1);
			int parentNode2 = parentNodeIdBelowFilter(childNode2);
			if(parentNode1 != parentNode2 || parentNode1 != i)
				__debugbreak();
		}
	}

	// Checks if the heap is valid
	void debugCheckHeap() {
		volatile int nodeId, parentNodeId, leftParentNodeId, rightParentNodeId;
		for(nodeId = 1; nodeId < m_last; nodeId++) {

			if(nodeId > filterLevelEnd) {
				parentNodeId = parentNodeIdBelowFilter(nodeId);
			} else if(nodeId < filterLevelStart + filterLevelSize) {
				parentNodeId = parentNodeIdAboveFilter(nodeId);
			} else {
				int filterLevelOffset = filterLevelUp(nodeId);
				if(filterConditionUp(nodeId, filterLevelOffset))
				{
					leftParentNodeId = parentNodeId1InFilter(nodeId);
					rightParentNodeId = rightParentNodeId2InFilter(nodeId, filterLevelOffset);
				} else {
					leftParentNodeId = leftParentNodeId2InFilter(nodeId, filterLevelOffset);
					rightParentNodeId = parentNodeId1InFilter(nodeId);
				}
				if(m_tmpError[m_nodes[correctNodeId(nodeId)]] < m_tmpError[m_nodes[correctNodeId(leftParentNodeId)]] || m_tmpError[m_nodes[correctNodeId(nodeId)]] < m_tmpError[m_nodes[correctNodeId(rightParentNodeId)]]) {
					volatile int vidc = m_nodes[correctNodeId(nodeId)];
					volatile int vidlp = m_nodes[correctNodeId(leftParentNodeId)];
					volatile int vidrp = m_nodes[correctNodeId(leftParentNodeId)];
					volatile float errtc = m_tmpError[vidc];
					volatile float errtlp = m_tmpError[vidlp];
					volatile float errtrp = m_tmpError[vidrp];
					__debugbreak();
				}
				continue;
			}

			if(m_tmpError[m_nodes[correctNodeId(nodeId)]] < m_tmpError[m_nodes[correctNodeId(parentNodeId)]]) {
				volatile int vidc = m_nodes[correctNodeId(nodeId)];
				volatile int vidp = m_nodes[correctNodeId(parentNodeId)];
				volatile float errtc = m_tmpError[vidc];
				volatile float errtp = m_tmpError[vidp];
				__debugbreak();
			}
		}
	}

protected:
	// This is an optimisation to prevent memory collisions (false sharing) in the first "m_countMemCollisionNodes" Nodes (and locks)
	const int m_countMemCollisionNodes = 191; // 128 + 63
	// Used to fetch the real nodeId with the "nodeId"
	inline int correctNodeId(int id) const {
		//return m_haltonIndices[id];
		if (id < m_countMemCollisionNodes) {
			return (16 * id);
		} 
		return (id + (m_countMemCollisionNodes) * 15);
	}

	inline int correctLockId(int id) const {
		if(id < m_countMemCollisionNodes) {
			return (16 * id);
		} 
		return (id + (m_countMemCollisionNodes) * 15);
	}

	int repairDown(int nodeId); // Iterative part of pop returns the next node or -1
	bool repairUp(int nodeId, int lookupId, int parentNodeId, int parentlookupId); // Assumes lock of nodeId and parent nodeId releases lock of nodId and keeps parent lock

	// The bubbleUp operation in the filter layers.
	bool repairUpFilter(int lookupId);

	int m_size;

	std::vector<CRITICAL_SECTION> m_locks; // Per Node Locks
	std::vector<std::atomic<int>> m_notRepairedUp; // Blocks the bubbleUp from child if the current node is not at the correct position
	std::vector<std::atomic<T>> m_tmpError; // An error cache used to maintain the structure
	std::vector<int> m_nodes;
	std::vector<std::atomic<int>> m_nodeLookup; // A lookup to check who is in which node
	std::atomic<int> m_last; // Index of the Node that is to be deleted on pop, the last element in the array;

	const int rootNodeIndex = 0; // Dont change this otherwise the other formulars have to be adjusted
	const int filteredElements = (filterLevelCount == 0 ? 0 : (2 << (filterLevelCount)) - 1); // 2^(filterLevelCount+1) - 1
	const int filterLevelSize = (filterLevelCount == 0 ? 0: (2 << (filterLevelCount-1)));	      // 2^(filterLevelCount) // or 0 if 0 filterLevels

	const int filterLevelStart = filteredElements - filterLevelSize;
	const int filterLevelEnd = filteredElements + (filterLevelSize * filterLevelCount) - 1;

	// Calculate the filter level offset if you go up or down respectively
	inline int filterLevelUp(int nodeId) const {
		int filterLevelOffset;
		filterLevelOffset = ((filterLevelEnd - nodeId) / filterLevelSize);
		//filterLevelOffset = ((nodeId - filterLevelStart) / filterLevelSize) - 1;
		filterLevelOffset = (1 << filterLevelOffset);
		return filterLevelOffset;
	}

	inline int filterLevelDown(int nodeId) const {
		int filterLevelOffset;
		filterLevelOffset = ((filterLevelEnd - nodeId) / filterLevelSize) - 1;
		//filterLevelOffset = ((nodeId - filterLevelStart) / filterLevelSize);
		filterLevelOffset = (1 << filterLevelOffset);
		return filterLevelOffset;
	}

	// Condition for the second parent/child node
	inline bool filterConditionUp(int nodeId, int filterLevelOffset) const {
		return(((nodeId - filterLevelStart) % (filterLevelSize)) < filterLevelOffset);
	}

	inline int filterConditionDown(int nodeId, int filterLevelOffset) const {
		return(((nodeId - filterLevelStart) % (filterLevelSize)) >= filterLevelSize - filterLevelOffset);
	}
	// Returns the corresponding index of left/right child/parent
	inline int leftChildNodeIdAboveFilter(int nodeId) const {
		return (nodeId * 2) + 1;
	}

	inline int rightChildNodeIdAboveFilter(int nodeId) const {
		return (nodeId * 2) + 2;
	}

	inline int parentNodeIdAboveFilter(int nodeId) const {
		return ((nodeId - 1) / 2);
	}

	// In Filter

	inline int childNodeId1InFilter(int nodeId) const {
		return (nodeId + filterLevelSize);
	}

	inline int parentNodeId1InFilter(int nodeId) const {
		return (nodeId - filterLevelSize);
	}

	inline int rightChildNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId + filterLevelSize + filterLevelOffset);
	}

	inline int leftChildNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId + filterLevelOffset);
	}

	inline int rightParentNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId - filterLevelOffset);
	}

	inline int leftParentNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId - filterLevelSize - filterLevelOffset);
	}

	// If you want the filter level connections to be the other way around you can use this instead (also have to use the commented filterLevelOffset)

	/*inline int rightChildNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId + filterLevelSize + filterLevelOffset);
	}

	inline int leftChildNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId + filterLevelSize - filterLevelOffset);
	}

	inline int rightParentNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId - filterLevelSize + filterLevelOffset);
	}

	inline int leftParentNodeId2InFilter(int nodeId, int filterLevelOffset) const {
		return (nodeId - filterLevelSize - filterLevelOffset);
	}*/

	// Below Filter

	inline int parentNodeIdBelowFilter(int nodeId) const {
		//return parentNodeIdAboveFilter(nodeId);
		return (((nodeId - 1  + (filterLevelSize * filterLevelCount)) / 2));
	}
	inline int leftChildNodeIdBelowFilter(int nodeId) const { 
		//return leftChildNodeIdAboveFilter(nodeId);
		return (nodeId * 2 + 1 - (filterLevelSize * filterLevelCount));
	}
	inline int rightChildNodeIdBelowFilter(int nodeId) const { 
		//return rightChildNodeIdAboveFilter(nodeId);
		return (nodeId * 2 + 2 - (filterLevelSize * filterLevelCount));
	}
};

template<typename T>
inline FPPQ<T>::FPPQ(int maxCapacity) {
	if(maxCapacity % 2)
		m_size = maxCapacity;
	else
		m_size = maxCapacity + 1; // We need the size to be odd
	m_locks = std::vector<CRITICAL_SECTION>(m_size + m_countMemCollisionNodes*15);
	m_notRepairedUp = std::vector<std::atomic<int>>(m_size);
	m_nodes = std::vector<int>(m_size + m_countMemCollisionNodes * 15,-1);
	m_tmpError = std::vector<std::atomic<T>>(m_size);
	m_nodeLookup = std::vector<std::atomic<int>>(m_size);
	int spincount1 = 0x01001000;
	for(int i = 0; i < m_size + m_countMemCollisionNodes * 15; i++) {
		if(!InitializeCriticalSectionAndSpinCount(&m_locks[i],spincount1))
			__debugbreak();
	}
	for(int i = 0; i < m_size; i++) {
		m_nodeLookup[i] = -1;
	}
	m_last = -1;
}

template<typename T>
void FPPQ<T>::fill(std::vector<std::pair<T, int>> errorValues, int countThreads) {
	std::sort(std::execution::par_unseq, errorValues.begin(), errorValues.end(), std::less<std::pair<T, int>>());
#pragma omp parallel for num_threads(countThreads)
	for(int i = 0; i < errorValues.size(); i++) {
		m_nodes[correctNodeId(i)] = errorValues[i].second;
		m_tmpError[errorValues[i].second] = errorValues[i].first;
		//m_nodes[i].setDebugError(errorValues[i].first);
		m_nodeLookup[errorValues[i].second] = i;

	}
	m_last = static_cast<int>(errorValues.size()) -1;
}

template<typename T>
int FPPQ<T>::pop() {
	int popIndex = 0;
	int popIncreaseCount = maxPopIndex;
	// The PopLevelBase Stuff is just to diverge the threads a little bit
	int popLevel = 0;
	int popLevelBase = 0;
	int popLevelIndex = 0;
	int lookupId = -1;
	while(lookupId == -1) {
		if(popIncreaseCount >= maxPopIndex) {
			popIndex = 0;
			popLevel = 0;
			popIncreaseCount = 1;
			popLevelBase = 0;
			popLevelIndex = 0;
		} else {
			if(popLevelIndex >= (1 << (popLevel)) - 1) {
				popLevel++;
				popIncreaseCount++;
				popLevelIndex = 0;
				popLevelBase = (1 << (popLevel)) - 1;
				popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
			} else {
				popLevelIndex++;
				popIncreaseCount++;
				popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
			}
			

		}
		// Lock a root node
		while(lookupId == -1) {
			while(!tryLock(popIndex))  // acquire lock
			{
				if(popIncreaseCount >= maxPopIndex) {
					popIndex = 0;
					popLevel = 0;
					popIncreaseCount = 1;
					popLevelBase = 0;
					popLevelIndex = 0;
				} else {
					if(popLevelIndex >= (1 << (popLevel)) - 1) {
						popLevel++;
						popIncreaseCount++;
						popLevelIndex = 0;
						popLevelBase = (1 << (popLevel)) - 1;
						popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
					} else {
						popLevelIndex++;
						popIncreaseCount++;
						popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
					}

				}
			}
			lookupId = m_nodes[correctNodeId(popIndex)];
			if(lookupId == -1) {
				unlock(popIndex);
				//__debugbreak();
				if(m_last < 0)
					return -1;
				popIndex = 0;
				popLevel = 0;
				popIncreaseCount = 1;
				popLevelBase = 0;
				popLevelIndex = 0;
			}
			// This should not happen
			else if(m_nodeLookup[lookupId] == -1) {
				unlock(popIndex);
				__debugbreak();
				popIncreaseCount++;
				lookupId = -1;
			}
		}

		{
			int lastIndex = m_last;
			int oldIndex = lastIndex;
			bool failed = false;
			if(lastIndex >= popIndex && lastIndex != -1) {
				while(!m_last.compare_exchange_weak(oldIndex, lastIndex - 1, std::memory_order_acquire)) {
					lastIndex = m_last;
					oldIndex = lastIndex;
					if(lastIndex < popIndex || lastIndex == -1) {
						failed = true;
						break;
					}
				}
			} else {
				failed = true;
			}
			if(failed) {
				unlock(popIndex);
				lookupId = -1;
			}
			else {
				m_nodeLookup[lookupId].store(-1);
				if(lastIndex > popIndex) {
					lock(lastIndex);
					int lastLookupId = m_nodes[correctNodeId(lastIndex)];
					while(lastLookupId == -1) {
						unlock(lastIndex);
						std::this_thread::yield();
						lock(lastIndex);
						lastLookupId = m_nodes[correctNodeId(lastIndex)];
					}

					m_nodeLookup[lastLookupId].store(popIndex);
					m_nodes[correctNodeId(popIndex)] = (lastLookupId);
					m_notRepairedUp[lastLookupId] += 1;
					m_nodes[correctNodeId(lastIndex)] = (-1);
					unlock(lastIndex);

					int nextNodeId = popIndex;
					while(nextNodeId != -1) {
						nextNodeId = repairDown(nextNodeId);
					}

					m_notRepairedUp[lastLookupId] -= 1;
				} else { 
					if(lastIndex == popIndex) {
						m_nodes[correctNodeId(lastIndex)] = (-1);
						unlock(lastIndex);
					} else {
						// This should be impossible we are popping an element that is behind our last Index which can happen (because of multi-Threading) but should be prevented with the atomic CAS
						__debugbreak();
						throw std::exception("Invalid pop: popIndex after lastIndex");
					}
				}
			}
		}
		popIncreaseCount++;
	}

	return lookupId;

}


template<typename T>
void FPPQ<T>::update(int lookupId, T errorValue, bool fromInsertion) {
	bool directionUp = false;

	if(!fromInsertion) {
		T oldError = m_tmpError[lookupId];
		if(oldError == errorValue) {
			return;
		} else {
			m_notRepairedUp[lookupId] += 1; // mark as closed for vertices traversing in upwards direction
			if(errorValue < oldError) { // We now have a smaller Error (We need to look upwards)
				directionUp = true;
				m_tmpError[lookupId] = errorValue; // TODO: this is not right could potentially lead to corrupted heap (does not happen in practice maybe this is ok after all?)
			}
		}
	} else {
		directionUp = true;
	}



	bool finished = false;
	int nodeId;
	if(directionUp) {
		int parentNodeId;
		int parentLookupId;
		while(!finished) {
			nodeId = m_nodeLookup[lookupId];

			if(nodeId <= rootNodeIndex) { // We are finished because node is one of the root nodes
				finished = true;
			} else {
				if(nodeId > filterLevelEnd) {
					parentNodeId = parentNodeIdBelowFilter(nodeId);
					if(filterLevelCount == 0) {
						if(nodeId == 0) {
							finished = true;
						}
					}
				} else if(nodeId < filterLevelStart + filterLevelSize) {
					if(nodeId <= 0) { // We are finished because node is the one with the smallest error or lookup got popped
						finished = true;
					} else {
						parentNodeId = parentNodeIdAboveFilter(nodeId);
					}
				} else { // We are in the Filter Levels
					if(repairUpFilter(lookupId)) {
						finished = true;
					}
					parentNodeId = parentNodeIdAboveFilter(nodeId);
				}
				parentLookupId = m_nodes[correctNodeId(parentNodeId)];
				while(!finished && (parentLookupId == -1 || m_notRepairedUp[parentLookupId])) // Critical for time saving check if parent is not up to date
				{
					std::this_thread::yield(); // saves time when we have too many threads
					nodeId = m_nodeLookup[lookupId];
					if(nodeId > filterLevelEnd) {
						parentNodeId = parentNodeIdBelowFilter(nodeId);
						if(filterLevelCount == 0) {
							if(nodeId == 0) {
								finished = true;
							}
						}
					} else if(nodeId < filterLevelStart + filterLevelSize) {
						if(nodeId <= 0) { // We are finished because node is the one with the smallest error or element got popped
							finished = true;
						} else {
							parentNodeId = parentNodeIdAboveFilter(nodeId);
						}
					} else {
						if(repairUpFilter(lookupId)) {
							finished = true;
						}
						parentNodeId = parentNodeIdAboveFilter(nodeId);
					}
					parentLookupId = m_nodes[correctNodeId(parentNodeId)];
				}
				if(!finished) {
					// Locking from parent to child is faster than the original implementation with a trylock after locking the child first.
					lock(parentNodeId);
					lock(nodeId);

					if(m_nodes[correctNodeId(nodeId)] != lookupId) { // Element got swapped to different node in the meantime
						unlock(nodeId);
						unlock(parentNodeId);
						// Try again
					} else {
						// upwards we need to compare the error with parents error
						bool stopInnerLoop = false;
						while(!(stopInnerLoop || finished)) { // inner loop (like the original implementation) fetch parents lock check if it can be swapped with the child node. if tryLock failed start again from the outer loop.
															  // Invariant: have lock of parentNodeId and nodeId 
							parentLookupId = m_nodes[correctNodeId(parentNodeId)];
							if(parentLookupId == -1 || m_notRepairedUp[parentLookupId]) { // If parent is currently closed
																// Try again
								unlock(nodeId);
								unlock(parentNodeId);
								stopInnerLoop = true;
								std::this_thread::yield();
							} else {
								// Dont clear parentNodeElement yet somebody may wants to close it.
								finished = repairUp(nodeId, lookupId, parentNodeId, parentLookupId);
								nodeId = parentNodeId;
								if(!finished) {
									if(nodeId > filterLevelEnd) {
										parentNodeId = parentNodeIdBelowFilter(nodeId);
										if(!tryLock(parentNodeId)) {
											unlock(nodeId);
											stopInnerLoop = true;
										}
									} else if(nodeId < filterLevelStart + filterLevelSize) {
										parentNodeId = parentNodeIdAboveFilter(nodeId);
										if(!tryLock(parentNodeId)) {
											unlock(nodeId);
											stopInnerLoop = true;
										}
									} else {
										unlock(nodeId);
										if(repairUpFilter(lookupId)) {
											finished = true;
										}
										stopInnerLoop = true;
									}
								} else {
									unlock(nodeId);
								}
							}
						}
					}
				}
			}
		}
	} else { // downwards
		while(!finished) {
			nodeId = m_nodeLookup[lookupId];
			if(nodeId < 0) { // We are finished because node got popped in theory this should not happen (at least in mesh decimation. could happen in other use cases)
				finished = true;
				//__debugbreak();
			} else {
				lock(nodeId);
				if(m_nodes[correctNodeId(nodeId)] != lookupId) { // Element got swapped to different node in the meantime
					unlock(nodeId);
					// Try again
				} else {
					m_tmpError[lookupId] = errorValue;
					int nextNodeId = nodeId;
					while(nextNodeId != -1)
						nextNodeId = repairDown(nextNodeId); // downwards like 
					finished = true;
				}
			}
		}
	}
	m_notRepairedUp[lookupId] -= 1; // mark as not closed anymore
}

template<typename T>
void FPPQ<T>::insert(int lookupId, T errorValue) {

	if(m_nodeLookup[lookupId] > 0) // Already inside queue
		__debugbreak();

	m_notRepairedUp[lookupId] += 1;
	m_tmpError[lookupId] = errorValue;

	int lastIndex = ++m_last;
	lock(lastIndex);

	while(m_nodes[correctNodeId(lastIndex)] != -1) {
		unlock(lastIndex);
		std::this_thread::yield();
		lock(lastIndex);
	}

	m_nodeLookup[lookupId].store(lastIndex);
	m_nodes[correctNodeId(lastIndex)] = lookupId;
	unlock(lastIndex);

	update(lookupId, errorValue, true);
}

template<typename T>
int FPPQ<T>::repairDown(int nodeId) // Assumes we have lock on nodeId
{
	if(nodeId < 0)
		__debugbreak();
	int nextNodeId = -1;
	int leftChildNodeId, rightChildNodeId;

	if(nodeId > filterLevelEnd - filterLevelSize) {
		leftChildNodeId = leftChildNodeIdBelowFilter(nodeId);
		rightChildNodeId = rightChildNodeIdBelowFilter(nodeId);
	} else if(nodeId < filterLevelStart) {
		leftChildNodeId = leftChildNodeIdAboveFilter(nodeId);
		rightChildNodeId = rightChildNodeIdAboveFilter(nodeId);
	} else {
		int filterLevelOffset = filterLevelDown(nodeId);
		if(filterConditionDown(nodeId, filterLevelOffset))
		{
			leftChildNodeId = leftChildNodeId2InFilter(nodeId, filterLevelOffset);
			rightChildNodeId = childNodeId1InFilter(nodeId);
		} else {
			leftChildNodeId = childNodeId1InFilter(nodeId);
			rightChildNodeId = rightChildNodeId2InFilter(nodeId, filterLevelOffset);
		}
	}

	if(leftChildNodeId <= m_last) { // We assume an odd number of total nodes. Could be dangerous
		lock(leftChildNodeId);
		lock(rightChildNodeId);
		{ 
			T valueLeft, valueRight;
			int lookupIdLeft = m_nodes[correctNodeId(leftChildNodeId)];
			int lookupIdRight = m_nodes[correctNodeId(rightChildNodeId)];

			assert(lookupIdLeft != lookupIdRight || lookupIdLeft == -1);

			if(lookupIdLeft == -1) {
				valueLeft = (std::numeric_limits<T>::max)();
			} else {
				valueLeft = m_tmpError[lookupIdLeft];
			}
			if(lookupIdRight == -1) {
				valueRight = (std::numeric_limits<T>::max)();
			} else {
				valueRight = m_tmpError[lookupIdRight];
			}
			int lookupId = m_nodes[correctNodeId(nodeId)];
			if(valueLeft >= valueRight) {
				if(valueRight < m_tmpError[lookupId]) {
					nextNodeId = rightChildNodeId;
					unlock(leftChildNodeId); // release the other child

					int lookupIdChild = lookupIdRight;

					m_nodeLookup[lookupId].store(nextNodeId); // swap
					m_nodes[correctNodeId(nextNodeId)] = (lookupId);
					m_nodeLookup[lookupIdChild].store(nodeId);
					m_nodes[correctNodeId(nodeId)] = (lookupIdChild);
				} else {
					// we have the correct constellation
					// release the childs
					unlock(leftChildNodeId);
					unlock(rightChildNodeId);

				}
			} else if(valueLeft < m_tmpError[lookupId]) {
				nextNodeId = leftChildNodeId;
				unlock(rightChildNodeId); // release the other child

				int lookupIdChild = lookupIdLeft;
				m_nodeLookup[lookupId].store(nextNodeId); // swap
				m_nodes[correctNodeId(nextNodeId)] = (lookupId);
				m_nodeLookup[lookupIdChild].store(nodeId);
				m_nodes[correctNodeId(nodeId)] = (lookupIdChild);
			} else {
				// we have the correct constellation
				unlock(leftChildNodeId);
				unlock(rightChildNodeId); // release the childs
			}
		} 
	}
	unlock(nodeId); // release Node because we are finished

	return nextNodeId; // Invariant We own lock to nextNodeId (if it exists) and have released the lock of nodeId
	
}




template<typename T>
bool FPPQ<T>::repairUp(int nodeId, int lookupId, int parentNodeId, int parentlookupId) // Assumes lock of nodeId, parent nodeId und parent lookupId releases lock of nodId and keeps parent node lock
{
	float valueParent = m_tmpError[parentlookupId];
	if(nodeId < 0)
		__debugbreak();
	if(m_tmpError[lookupId] < valueParent) { // We are now smaller than the parent
		// Node
		m_nodeLookup[lookupId].store(parentNodeId); // swap
		m_nodes[correctNodeId(parentNodeId)] = (lookupId);
		// Parent
		m_nodeLookup[parentlookupId].store(nodeId);
		m_nodes[correctNodeId(nodeId)] = (parentlookupId);
		// Locks
		unlock(nodeId);
		if(nodeId <= ((rootNodeIndex + 1) * 3) - 1)
			return true; // The node is now the root so we are finished
		return false; // We have to look at our parent
	} else {
		// we have the correct constellation
		unlock(nodeId);
		return true;
	}
}

template<typename T>
bool FPPQ<T>::repairUpFilter(int lookupId) {
	bool finished = false;
	int nodeId, leftParentNodeId, rightParentNodeId, leftParentLookupId, rightParentLookupId;
	int nextNodeId = -1;
	while(!finished) {
		nodeId = m_nodeLookup[lookupId];
		if(nodeId < filterLevelStart + filterLevelSize) { // We are over the filter Levels
			return false;
		} else {
			if(nodeId > filterLevelEnd) // an other thread updated us in downwards direction
			{
				return false;
			}

			int filterLevelOffset = filterLevelUp(nodeId);
			if(filterConditionUp(nodeId, filterLevelOffset))
			{
				leftParentNodeId = parentNodeId1InFilter(nodeId);
				rightParentNodeId = rightParentNodeId2InFilter(nodeId, filterLevelOffset);
			} else {
				leftParentNodeId = leftParentNodeId2InFilter(nodeId, filterLevelOffset);
				rightParentNodeId = parentNodeId1InFilter(nodeId);
			}
			leftParentLookupId = m_nodes[correctNodeId(leftParentNodeId)];
			rightParentLookupId = m_nodes[correctNodeId(rightParentNodeId)];
			while(!finished && (leftParentLookupId == -1 || rightParentLookupId == -1 || m_notRepairedUp[leftParentLookupId] || m_notRepairedUp[rightParentLookupId])) // Critical for time saving
			{
				std::this_thread::yield(); // saves time when we have too many threads
				nodeId = m_nodeLookup[lookupId];
				if(nodeId < filterLevelStart + filterLevelSize) { // we are now over the filter levels
					return false;
				} else {
					if(nodeId > filterLevelEnd) // an other thread updated us in downwards direction
					{
						return false;
					}
					int filterLevelOffset = filterLevelUp(nodeId);
					if(filterConditionUp(nodeId, filterLevelOffset)) 
					{
						leftParentNodeId = parentNodeId1InFilter(nodeId);
						rightParentNodeId = rightParentNodeId2InFilter(nodeId, filterLevelOffset);;
					} else {
						leftParentNodeId = leftParentNodeId2InFilter(nodeId, filterLevelOffset);
						rightParentNodeId = parentNodeId1InFilter(nodeId);
					}
				}
				leftParentLookupId = m_nodes[correctNodeId(leftParentNodeId)];
				rightParentLookupId = m_nodes[correctNodeId(rightParentNodeId)];
			}
			if(!finished) {
				lock(leftParentNodeId);
				bool aquired = true;
				int counter = 0;
				while(!tryLock(rightParentNodeId)) {  // acquire lock right child
					counter++;
					if(counter > 16) {
						aquired = false;
						break;
					}
				}
				if(!aquired) {
					unlock(leftParentNodeId);
				} else {
					float valueLeft, valueRight;
					leftParentLookupId = m_nodes[correctNodeId(leftParentNodeId)];
					rightParentLookupId = m_nodes[correctNodeId(rightParentNodeId)];

					assert(leftParentLookupId != rightParentLookupId);

					if(leftParentLookupId == -1 || rightParentLookupId == -1 || m_notRepairedUp[leftParentLookupId] || m_notRepairedUp[rightParentLookupId]) {
						unlock(leftParentNodeId);
						unlock(rightParentNodeId);
					} else {
						lock(nodeId);// acquire lock self
						if(m_nodes[correctNodeId(nodeId)] != lookupId) // Element got swapped to different position
						{
							unlock(nodeId);
							unlock(leftParentNodeId);
							unlock(rightParentNodeId);
						} else {
							valueLeft = m_tmpError[leftParentLookupId];
							valueRight = m_tmpError[rightParentLookupId];
							int lookupIdDebug = m_nodes[correctNodeId(nodeId)];
							if(valueLeft <= valueRight) {
								if(valueRight > m_tmpError[lookupId]) { // If the bigger parent is bigger than myself
									nextNodeId = rightParentNodeId;
									unlock(leftParentNodeId); // release the other child

									int lookupIdParent = rightParentLookupId;

									m_nodeLookup[lookupId].store(nextNodeId); // swap
									m_nodes[correctNodeId(nextNodeId)] = (lookupId);
									m_nodeLookup[lookupIdParent].store(nodeId);
									m_nodes[correctNodeId(nodeId)] = (lookupIdParent);
									unlock(rightParentNodeId); // release the childs
									unlock(nodeId); // release the childs
								} else {
									// we have the correct constellation
									unlock(nodeId);
									unlock(leftParentNodeId);
									unlock(rightParentNodeId);
									finished = true;
								}
							} else if(valueLeft > m_tmpError[lookupId]) {
								nextNodeId = leftParentNodeId;
								unlock(rightParentNodeId); // release the other child

								int lookupIdParent = leftParentLookupId;

								m_nodeLookup[lookupId].store(nextNodeId); // swap
								m_nodes[correctNodeId(nextNodeId)] = (lookupId);
								m_nodeLookup[lookupIdParent].store(nodeId);
								m_nodes[correctNodeId(nodeId)] = (lookupIdParent);
								unlock(leftParentNodeId); // release the childs
								unlock(nodeId); // release the childs
							} else {
								// we have the correct constellation
								unlock(nodeId);
								unlock(leftParentNodeId);
								unlock(rightParentNodeId); // release the childs
								finished = true;
							}
						}
					}

				}
			}
		}
	}
	return true; // Invariant No locks
}