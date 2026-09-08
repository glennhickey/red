/*
 * ChromListMaker.h
 *
 *  Created on: Mar 13, 2014
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef CHROMLISTMAKER_H_
#define CHROMLISTMAKER_H_

#include <string>
#include <vector>

#include "Chromosome.h"
#include "ChromosomeOneDigit.h"

#include "../utility/Util.h"

using namespace std;
using namespace utility;

namespace nonltr {

class ChromListMaker {
private:
	vector<Chromosome *> * chromList;
	string seqFile;

	// State of the one-sequence-at-a-time reader.  Named apart from the local
	// "in" of the two list-building readers, which are independent of it.
	ifstream streamIn;
	bool isStreamOpen;
	string pendingHeader;

public:
	ChromListMaker(string);
	virtual ~ChromListMaker();
	const vector<Chromosome *> * makeChromList();
	const vector<Chromosome *> * makeChromOneDigitList();

	/**
	 * Read the next record of the file, one sequence at a time.
	 *
	 * The two methods above hold every sequence of the file in memory at once,
	 * which for a single-file genome means the whole genome -- and the scan
	 * stage builds two such lists, one encoded and one not, so it carries two
	 * complete copies.  Reading a record at a time bounds the cost by the
	 * longest sequence instead.
	 *
	 * Returns false at end of file.  hadSequence reports whether the record had
	 * any sequence lines at all; Chromosome treats a record with none as an
	 * error, and passing the flag on keeps that behaviour.
	 */
	bool nextSequence(string & header, string & seq, bool & hadSequence);

	/* Build one chromosome from a record returned by nextSequence. */
	static Chromosome * makeChrom(string & header, string & seq,
			bool hadSequence);
	static ChromosomeOneDigit * makeChromOneDigit(string & header, string & seq,
			bool hadSequence);
};

} /* namespace nonltr */
#endif /* CHROMLISTMAKER_H_ */
