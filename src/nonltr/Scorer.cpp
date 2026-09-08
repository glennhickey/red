/*
 * Scorer.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */
#include "Scorer.h"

Scorer::Scorer(ChromosomeOneDigit * chromIn,
		ITableView<unsigned long, int> * const table) {
	chrom = chromIn;
	kmerTable = table;
	scores = new vector<int>(chrom->getBase()->size(), 0);
	k = kmerTable->getK();
	max = -1;
	score();
	checkNotEmpty();
}

Scorer::~Scorer() {
	scores->clear();
	delete scores;
}

/**
 * This method scores each nucleotide in the chromosome.
 * The nucleotides represented by 'N' are assigned zero.
 */
void Scorer::score() {
	const vector<vector<int> *> * segment = chrom->getSegment();
	const char * segBases = chrom->getBase()->c_str();

	for (int s = 0; s < segment->size(); s++) {
		int start = segment->at(s)->at(0);
		int end = segment->at(s)->at(1);
		kmerTable->wholesaleValueOf(segBases, start, end - k + 1, scores,
				start);

		// Handle the last word from end - k + 2 till the end, inclusive.
		int * sc = scores->data();
		for (int i = end - k + 2; i <= end; i++) {
			sc[i] = sc[i - 1];
		}
	}
}

/**
 * The logarithm of every score that fits in the cache below, tabulated.
 *
 * takeLog runs over the whole genome four times in a -gnm run -- once while
 * training and once per strand while scanning -- and a libm log() per base is
 * a large share of that.  Scores are small counts, so almost all of them fall
 * inside the table; anything above it still goes through log().  The table
 * entries are produced by the identical expression, so the tabulated answers
 * are bit-for-bit what the call would have returned.
 *
 * The table is kept across calls because a fragmented assembly can hold
 * hundreds of thousands of sequences, and rebuilding per sequence would cost
 * more than it saves.
 */
namespace {
const int LOG_CACHE_SIZE = 1 << 16;

struct LogCache {
	double logBase;
	bool isBuilt;
	vector<int> value;

	LogCache() : logBase(0.0), isBuilt(false) {
	}

	const int * get(double logBaseIn) {
		if (!isBuilt || logBase != logBaseIn) {
			value.resize(LOG_CACHE_SIZE);
			// Entry 0 is never read: takeLog leaves a score of zero alone.
			value[0] = 0;
			for (int i = 1; i < LOG_CACHE_SIZE; i++) {
				value[i] = (int) ceil(log((double) i) / logBaseIn);
			}
			logBase = logBaseIn;
			isBuilt = true;
		}
		return value.data();
	}
};

// thread_local so that this stays correct if the per-chromosome work is ever
// run on more than one thread; the table is read once per takeLog call, not
// per base, so the TLS access costs nothing measurable.
thread_local LogCache logCache;
}

/**
 * This method takes the logarithm of the scores according to the base.
 * If the score equals zero, it is left the same.
 */
void Scorer::takeLog(double base) {
	// Handle the case where base is one
	bool isOne = false;
	if (fabs(base - 1.0) < std::numeric_limits<double>::epsilon()) {
		isOne = true;
	}
	double logBase = isOne ? log(1.5) : log(base);

	const int * const logOf = logCache.get(logBase);
	const int lowest = isOne ? 2 : 1;

	const vector<vector<int> *> * segment = chrom->getSegment();
	int * const sc = scores->data();
	for (int s = 0; s < segment->size(); s++) {
		int start = segment->at(s)->at(0);
		int end = segment->at(s)->at(1);
		for (int h = start; h <= end; h++) {
			int score = sc[h];

			// score == 0 is left alone, and so is score == 1 when the base was
			// adjusted up from one, exactly as before.
			if (score >= lowest) {
				sc[h] = (score < LOG_CACHE_SIZE) ?
						logOf[score] : (int) ceil(log((double) score) / logBase);
			}
		}
	}
}

int Scorer::getK() {
	return k;
}

vector<int>* Scorer::getScores() {
	return scores;
}

void Scorer::printScores(string outputFile, bool canAppend) {
	ofstream outScores;
	if (canAppend) {
		outScores.open(outputFile.c_str(), ios::out | ios::app);
	} else {
		outScores.open(outputFile.c_str(), ios::out);
	}

	Util::checkStream(outScores, outputFile, "open");

	int step = 50;
	outScores << chrom->getHeader() << endl;
	int len = scores->size();
	for (int i = 0; i < len; i = i + step) {
		int e = (i + step - 1 > len - 1) ? len - 1 : i + step - 1;
		for (int k = i; k <= e; k++) {
			outScores << scores->at(k) << " ";
		}
		outScores << endl;
	}
	outScores << endl;

	outScores.flush();
	Util::checkStream(outScores, outputFile, "write to");
	outScores.close();
	Util::checkStream(outScores, outputFile, "close");
}

int Scorer::countLessOrEqual(int thr) {
	int count = 0;
	const vector<vector<int> *> * segment = chrom->getSegment();
	const int * const sc = scores->data();
	for (int s = 0; s < segment->size(); s++) {
		int start = segment->at(s)->at(0);
		int end = segment->at(s)->at(1);
		for (int h = start; h <= end; h++) {
			if (sc[h] <= thr) {
				count++;
			}
		}
	}
	return count;
}

/**
 * The constructor used to sweep every scored base to find the maximum, but
 * nothing ever read it -- getMax() has no caller inside Red.  That sweep was a
 * whole extra pass over the genome per Scorer, and five Scorers are built over
 * the course of a -gnm run.  All the pass could actually detect was a
 * chromosome with nothing scored in it, which is what this checks instead; the
 * maximum itself is computed on demand.
 */
void Scorer::checkNotEmpty() {
	const vector<vector<int> *> * segmentList = chrom->getSegment();
	int segmentCount = segmentList->size();
	for (int jj = 0; jj < segmentCount; jj++) {
		vector<int> * segment = segmentList->at(jj);
		if (segment->at(1) >= segment->at(0)) {
			return;
		}
	}

	string msg("Error occurred while finding the maximum score.");
	throw InvalidStateException(msg);
}

int Scorer::getMax() {
	if (max == -1) {
		const vector<vector<int> *> * segmentList = chrom->getSegment();
		int segmentCount = segmentList->size();
		const int * const sc = scores->data();
		for (int jj = 0; jj < segmentCount; jj++) {
			vector<int> * segment = segmentList->at(jj);
			int start = segment->at(0);
			int end = segment->at(1);
			for (int ss = start; ss <= end; ss++) {
				if (sc[ss] > max) {
					max = sc[ss];
				}
			}
		}
	}
	return max;
}
