/*
 * ChromListMaker.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Hani Zakaira Girgis
 */

#include "ChromListMaker.h"

namespace nonltr {

ChromListMaker::ChromListMaker(string seqFileIn) {
	seqFile = seqFileIn;
	chromList = new vector<Chromosome *>();
	isStreamOpen = false;
}

ChromListMaker::~ChromListMaker() {
	Util::deleteInVector(chromList);
	delete chromList;
}

const vector<Chromosome *> * ChromListMaker::makeChromList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	Chromosome * chrom = nullptr;

	while (in.good()) {
		string line;
		getline(in, line);
		if (line.empty()) {
			continue;  // Skip empty lines
		}
		if (line[0] == '>') {
			if (!isFirst) {
				chrom->finalize();
				chromList->push_back(chrom);
			} else {
				isFirst = false;
			}

			chrom = new Chromosome();
			chrom->setHeader(line);
		} else if (chrom != nullptr) {
			chrom->appendToSequence(line);
		}
	}
	if (chrom != nullptr) {
		chrom->finalize();
		chromList->push_back(chrom);
	}
	in.close();

	return chromList;
}

const vector<Chromosome *> * ChromListMaker::makeChromOneDigitList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	ChromosomeOneDigit * chrom = nullptr;

	while (in.good()) {
		string line;
		getline(in, line);
		if (line.empty()) {
			continue;  // Skip empty lines
		}
		if (line[0] == '>') {
			if (!isFirst) {
				chrom->finalize();
				chromList->push_back(chrom);
			} else {
				isFirst = false;
			}

			chrom = new ChromosomeOneDigit();
			chrom->setHeader(line);
		} else if (chrom != nullptr) {
			chrom->appendToSequence(line);
		}
	}

	if (chrom != nullptr) {
		chrom->finalize();
		chromList->push_back(chrom);
	}
	in.close();

	return chromList;
}

bool ChromListMaker::nextSequence(string & header, string & seq,
		bool & hadSequence) {
	if (!isStreamOpen) {
		streamIn.open(seqFile.c_str());
		isStreamOpen = true;

		// Anything before the first header is ignored, as it is by the
		// list-building readers above.
		string line;
		while (getline(streamIn, line)) {
			if (line.empty()) {
				continue;
			}
			if (line[0] == '>') {
				pendingHeader = line;
				break;
			}
		}
	}

	if (pendingHeader.empty()) {
		return false;
	}

	header = pendingHeader;
	pendingHeader.clear();

	// clear() rather than a fresh string, so the buffer earned by the longest
	// sequence so far is reused instead of being reallocated per record.
	seq.clear();
	hadSequence = false;

	string line;
	while (getline(streamIn, line)) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			pendingHeader = line;
			break;
		}
		seq.append(line);
		hadSequence = true;
	}

	return true;
}

Chromosome * ChromListMaker::makeChrom(string & header, string & seq,
		bool hadSequence) {
	Chromosome * chrom = new Chromosome();
	chrom->setHeader(header);
	if (hadSequence) {
		chrom->appendToSequence(seq);
	}
	chrom->finalize();
	return chrom;
}

ChromosomeOneDigit * ChromListMaker::makeChromOneDigit(string & header,
		string & seq, bool hadSequence) {
	ChromosomeOneDigit * chrom = new ChromosomeOneDigit();
	chrom->setHeader(header);
	if (hadSequence) {
		chrom->appendToSequence(seq);
	}
	chrom->finalize();
	return chrom;
}

}
/* namespace nonltr */
