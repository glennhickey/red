/*
 * TableBuilder.h
 *
 *  Created on: Jul 31, 2012
 *      Author: Hani Zakaria Girgis, PhD - NCBI/NLM/NIH
 */

#ifndef TABLEBUILDER_H_
#define TABLEBUILDER_H_

#include "KmerHashTable.h"
#include "EnrichmentMarkovView.h"
#include "ChromosomeOneDigit.h"
#include "ChromListMaker.h"
#include "IChromosome.h"

#include "../utility/Util.h"
#include "../exception/InvalidStateException.h"

#include <iostream>

using namespace std;
using namespace nonltr;
using namespace utility;
using namespace exception;

namespace nonltr {
class TableBuilder {
private:
	/**
	 * k-mer table
	 */
	KmerHashTable<unsigned long,int> * kmerTable;
	int maxValue;

	/**
	 * Directory including the FASTA files comprising the genome.
	 * These files must have the
	 */
	string genomeDir;

	/**
	 * The size of the motif
	 */
	int k;

	/**
	 * The total length of the whole genome
	 */
	long genomeLength;

	/**
	 * Count cap (-cap).  A percentile of the counts of the distinct k-mers in the
	 * genome, both strands pooled; every k-mer counted above the count at that
	 * percentile is flagged in highKmers, one bit per table entry, and the
	 * scanner masks every position such a k-mer covers.  capPct <= 0 turns it off.
	 */
	double capPct;
	long capThreshold;
	vector<unsigned long> highKmers;

	/**
	 * Methods
	 */
	void buildTable();
	void updateTable(ChromosomeOneDigit *);
	void findHighKmers();

public:
	TableBuilder(string, int, int, int, double capPct = 0.0);
	virtual ~TableBuilder();
	KmerHashTable<unsigned long,int> * const getKmerTable();
	void printTable();
	long getGenomeLength();
	int getMaxValue();
	bool hasCap() const;
	long getCapThreshold() const;
	const vector<unsigned long> & getHighKmers() const;
};
}

#endif /* TABLEBUILDER_H_ */
