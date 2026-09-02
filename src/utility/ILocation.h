/*
 * ILocation.h
 *
 *  Created on: Dec 20, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef ILOCATION_H_
#define ILOCATION_H_

#include <string>

using namespace std;

namespace utility {

class ILocation {
public:
	/*
	 * Required, not merely good practice: Location and EmptyLocation objects are
	 * routinely deleted through an ILocation*, for example by
	 * Util::deleteInVector<ILocation> and in LocationList and Scanner.  Without a
	 * virtual destructor here that is undefined behaviour, and it is not benign --
	 * the compiler emits the sized operator delete(void*, sizeof(ILocation)), which
	 * is smaller than the objects actually allocated.  glibc's allocator ignores
	 * that size and recovers the true one from its chunk header, so the bug hides;
	 * jemalloc trusts it, frees into the wrong size class and corrupts the heap,
	 * which resurfaces far away as scrambled Location bounds and a spurious
	 * "This list is not sorted" failure.
	 */
	virtual ~ILocation() {}

	virtual int getEnd() const = 0;
	virtual int getStart() const = 0;
	virtual void setEnd(int) = 0;
	virtual void setStart(int) = 0;
	virtual int getLength() = 0;
	virtual string toString() = 0;
};

}

#endif /* ILOCATION_H_ */
