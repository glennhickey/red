/*
 * TableBuilder.cpp
 *
 *  Created on: Jul 31, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "TableBuilder.h"

#include <cmath>
#include <cstdint>

TableBuilder::TableBuilder(string dir, int motifSize, int order, int minObs,
		double capPctIn) {
	genomeDir = dir;
	k = motifSize;
	genomeLength = 0;
	maxValue = 0;  // Initialize to prevent undefined behavior
	capPct = capPctIn;
	capThreshold = -1;
	// kmerTable = new KmerHashTable(k);
	// kmerTable = new EnrichmentView(k);

	// Whenever you change the template, modify line 50 and 70 and the header file line 35
	kmerTable = new EnrichmentMarkovView<unsigned long, int>(k, order, minObs);

	buildTable();
}

TableBuilder::~TableBuilder() {
	delete kmerTable;
}

void TableBuilder::buildTable() {
	vector<string> * fileList = new vector<string>();
	Util::readChromList(genomeDir, fileList, "fa");

	// One sequence in memory at a time rather than the whole file at once.
	string header;
	string seq;
	bool hadSequence = false;

	for (int i = 0; i < fileList->size(); i++) {
		cout << "Counting k-mers in " << fileList->at(i) << " ..." << endl;
		ChromListMaker maker(fileList->at(i));

		while (maker.nextSequence(header, seq, hadSequence)) {
			ChromosomeOneDigit * chrom = ChromListMaker::makeChromOneDigit(
					header, seq, hadSequence);
			genomeLength += chrom->getEffectiveSize();
			updateTable(chrom);
			delete chrom;
		}
	}
	// Check if overflow has occurred
	kmerTable->checkOverflow();

	// The raw counts are only in the table until processTable replaces them with
	// adjusted ones, so the count cap has to be worked out here.
	if (capPct > 0) {
		findHighKmers();
	}

	// View
	// EnrichmentView * view = dynamic_cast<EnrichmentView *>(kmerTable);
	EnrichmentMarkovView<unsigned long, int> * view =
			dynamic_cast<EnrichmentMarkovView<unsigned long, int> *>(kmerTable);

	if (view) {
		view->generateProbapilities();
		view->processTable();
		maxValue = view->getMaxValue();
	} else {
		throw InvalidStateException(string("Dynamic cast failed."));
	}
	cout << "Enrichment view is ready." << endl;

	fileList->clear();
	delete fileList;

	/* If you would like to see the contents of the table.*/
	// kmerTable-> printTable();
}

void TableBuilder::updateTable(ChromosomeOneDigit * chrom) {
	// EnrichmentView * view = dynamic_cast<EnrichmentView *>(kmerTable);
	EnrichmentMarkovView<unsigned long, int> * view =
			dynamic_cast<EnrichmentMarkovView<unsigned long, int> *>(kmerTable);

	const vector<vector<int> *> * segment = chrom->getSegment();
	const char * segBases = chrom->getBase()->c_str();

	for (int s = 0; s < segment->size(); s++) {
		int start = segment->at(s)->at(0);
		int end = segment->at(s)->at(1);
		// cerr << "The segment length is: " << (end-start+1) << endl;

		// Fast, but require some memory proportional to the segment length.
		kmerTable->wholesaleIncrement(segBases, start, end - k + 1);
		if (view) {
			view->count(segBases, start, end);
		} else {
			throw InvalidStateException(string("Dynamic cast failed."));
		}

		// Slow, but memory efficient
		/*
		 vector<int> hashList = vector<int>();
		 kmerTable->hash(segBases, start, end - k + 1, &hashList);

		 for (int i = start; i <= end - k + 1; i++) {
		 kmerTable->increment(segBases, i);
		 }
		 */
	}
}

KmerHashTable<unsigned long, int> * const TableBuilder::getKmerTable() {
	return kmerTable;
}

long TableBuilder::getGenomeLength() {
	if (genomeLength < 0) {
		string msg("The length of the genome cannot be negative.");
		throw InvalidStateException(msg);
	}

	return genomeLength;
}

int TableBuilder::getMaxValue() {
	return maxValue;
}

/**
 * The reverse complement of a k-mer key.  Keys hold the first base in the most
 * significant of k two-bit digits, A=0 C=1 G=2 T=3, so complementing is 3 - digit,
 * an xor, and reversing is reversing the order of the two-bit digits.
 */
static inline unsigned long reverseComplementKey(unsigned long key, int k) {
	uint32_t x = (uint32_t) key ^ (uint32_t) ((1UL << (2 * k)) - 1);
	x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
	x = ((x >> 4) & 0x0F0F0F0Fu) | ((x & 0x0F0F0F0Fu) << 4);
	x = ((x >> 8) & 0x00FF00FFu) | ((x & 0x00FF00FFu) << 8);
	x = (x >> 16) | (x << 16);
	return (unsigned long) (x >> (32 - 2 * k));
}

/**
 * Flag the k-mers whose count is above the capPct percentile of the counts of the
 * distinct k-mers in the genome.
 *
 * The table counts the forward strand only, so a k-mer's count on both strands is
 * its own entry plus its reverse complement's, and a k-mer and its reverse
 * complement are one distinct k-mer.  This is how WindowMasker counts its units,
 * and it too takes its thresholds as percentiles of that distribution, which is
 * what lets one setting serve genomes of any size.  The masking itself is not
 * WindowMasker's: WindowMasker clamps unit counts at t_high (by default the 99.8
 * percentile, so -cap 99.8 lands on the same count) and masks windows of units
 * whose average passes separate trigger and extension thresholds, where this
 * masks every k-mer above the cap outright.
 */
void TableBuilder::findHighKmers() {
	const int * values = kmerTable->getValues();
	const unsigned long size = kmerTable->getMaxTableSize();

	// Counts above histMax share the top bin.  The threshold only has to be exact
	// when it falls below histMax, and a percentile anywhere near the top of the
	// distribution lands far below it.
	const long histMax = 1L << 20;
	vector<long> hist(histMax + 1, 0);
	long distinct = 0;
	for (unsigned long x = 0; x < size; x++) {
		unsigned long y = reverseComplementKey(x, k);
		if (y < x) {
			continue;
		}
		long c = (long) values[x] + (y == x ? 0 : (long) values[y]);
		if (c > 0) {
			hist[c < histMax ? c : histMax]++;
			distinct++;
		}
	}

	long target = (long) ceil(capPct / 100.0 * distinct);
	long cumulative = 0;
	capThreshold = histMax;
	for (long c = 0; c <= histMax; c++) {
		cumulative += hist[c];
		if (cumulative >= target) {
			capThreshold = c;
			break;
		}
	}

	highKmers.assign(size / 64 + 1, 0);
	long flagged = 0;
	for (unsigned long x = 0; x < size; x++) {
		unsigned long y = reverseComplementKey(x, k);
		if (y < x) {
			continue;
		}
		long c = (long) values[x] + (y == x ? 0 : (long) values[y]);
		if (c > capThreshold) {
			highKmers[x >> 6] |= 1UL << (x & 63);
			highKmers[y >> 6] |= 1UL << (y & 63);
			flagged++;
		}
	}

	cout << "Count cap: the " << capPct << " percentile of " << distinct;
	cout << " distinct k-mers is a count of " << capThreshold << "; ";
	cout << flagged << " k-mers counted above it will be masked." << endl;
}

bool TableBuilder::hasCap() const {
	return capPct > 0;
}

long TableBuilder::getCapThreshold() const {
	return capThreshold;
}

const vector<unsigned long> & TableBuilder::getHighKmers() const {
	return highKmers;
}
