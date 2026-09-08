/*
 * EnrichmentMarkovView.cpp
 *
 *  Created on: Apr 17, 2013
 *      Author: Hani Zakaria Girgis, PhD
 */

namespace nonltr {

/**
 * The Markov order. It start at 0.
 */
template<class I, class V>
EnrichmentMarkovView<I, V>::EnrichmentMarkovView(int k, int order, int m) :
		minObs(m), factor(10000.00), KmerHashTable<I, V>(k) {
	initialize(order);
}

template<class I, class V>
EnrichmentMarkovView<I, V>::EnrichmentMarkovView(int k, V initValue, int order,
		int m) :
		minObs(m), factor(10000.00), KmerHashTable<I, V>(k, initValue) {
	initialize(order);
}

template<class I, class V>
void EnrichmentMarkovView<I, V>::initialize(int order) {
	// Test start
	// cout << "Testing: " << minObs << endl;
	// Test end

	o = order;
	if (o < 0) {
		string msg("The Markov order must be non-negative integer. ");
		msg.append("The invalid input is: ");
		msg.append(Util::int2string(o));
		msg.append(".");
		throw InvalidInputException(msg);
	}

	if (o >= KmerHashTable<I, V>::k) {
		string msg("The Markov order cannot be >= k (k-mer).");
		throw InvalidInputException(msg);
	}

	l = 0;
	modelList = new vector<KmerHashTable<int, int> *>();

	for (int i = 1; i <= o + 1; i++) {
		modelList->push_back(new KmerHashTable<int, int>(i));
	}
}

template<class I, class V>
EnrichmentMarkovView<I, V>::~EnrichmentMarkovView() {
	Util::deleteInVector(modelList);
	delete modelList;
}

/**
 * This method count words of size 1 to order+1 in the input sequence.
 * In other words, it updates the background tables. In addition, it
 * updates the length of the genome.
 *
 * sequence: is the input sequence.
 * start: the start index - inclosing.
 * end: the end index - inclosing.
 */
template<class I, class V>
void EnrichmentMarkovView<I, V>::count(const char * sequence, int start,
		int end) {

	// Multiple by 2 if scanning the forward strand and its reverse complement
	// l = l + (2 * (end - start + 1));
	l = l + (end - start + 1);

	int modelNumber = modelList->size();
	for (int i = 0; i < modelNumber; i++) {
		KmerHashTable<int, int> * t = modelList->at(i);
		t->wholesaleIncrement(sequence, start, end - i);
	}
}

/**
 * Normalize the count of words in each model.
 * Values stored in these models are multiplied by "factor."
 */
template<class I, class V>
void EnrichmentMarkovView<I, V>::generateProbapilities() {
	int modelNumber = modelList->size();

	for (int m = 0; m < modelNumber; m++) {
		KmerHashTable<int, int> * t = modelList->at(m);
		int tSize = t->getMaxTableSize();

		for (int i = 0; i < tSize; i += 4) {
			double sum = 0.0;

			for (int j = i; j < i + 4; j++) {
				sum += t->valueOf(j);
			}

			for (int j = i; j < i + 4; j++) {
				t->insert(j, round(factor * ((double) t->valueOf(j) / sum)));
			}
		}
	}
}

/**
 * Convert the raw k-mer counts into enrichment values.
 *
 * This visits every one of the 4^k table entries -- just over a billion at
 * k=15 -- so anything done per entry matters.  The original carried the key
 * around as a quaternary string, incremented it digit by digit, and asked the
 * background models to re-hash that string; for the highest-order model it
 * also allocated two vectors per key to collect the window values.
 *
 * The table index *is* the key: entry y holds the k-mer whose digits are the
 * base-4 digits of y, most significant first.  So every sub-word the models
 * are asked about is a shift and a mask of y, and no string, no hashing, and
 * no allocation is needed.  The arithmetic and the order of the multiplications
 * are unchanged.
 */
template<class I, class V>
void EnrichmentMarkovView<I, V>::processTable() {
	const int kLen = KmerHashTable<I, V>::k;
	const I tableSize = KmerHashTable<I, V>::maxTableSize;
	const int modelNumber = modelList->size();

	// The model of the highest order; its keys are o+1 bases long.
	KmerHashTable<int, int> * const oTable = modelList->at(modelNumber - 1);
	const int resultsSize = kLen - o - 1;
	const int wordMask = (1 << (2 * (o + 1))) - 1;

	double lowerP = 1.0;
	double upperP = 1.0;

	for (I y = 0; y < tableSize; y++) {
		if (y % 10000000 == 0) {
			cout << "Processing " << y << " keys out of "
					<< KmerHashTable<I, V>::maxTableSize;
			cout << endl;
		}

		// Calculate the expected number of occurrences.
		//
		// Both probabilities depend only on the first k-1 digits, so they are
		// shared by the four keys that differ in the last one.
		if (y % 4 == 0) {
			// a. Calculate probability from lower order models.
			lowerP = 1.0;
			for (int m = 0; m < modelNumber - 1; m++) {
				// The first m+1 digits of y.
				int prefix = (int) (y >> (2 * (kLen - m - 1)));
				lowerP *= (((double) modelList->at(m)->valueOf(prefix))
						/ factor);
			}

			// b. Calculate probability based on the specified order: the
			// windows of o+1 bases starting at 0 .. resultsSize-1.
			upperP = 1.0;
			for (int i = 0; i < resultsSize; i++) {
				int word = (int) ((y >> (2 * (kLen - i - o - 1))) & wordMask);
				upperP *= (((double) oTable->valueOf(word)) / factor);
			}
		}

		// The last window of o+1 bases, i.e. the one starting at resultsSize.
		const int lastWord = (int) (y & wordMask);

		// The expected number of occurances
		double exp = l * lowerP * upperP
				* (((double) oTable->valueOf(lastWord)) / factor);

		// Calculate the enrichment value.
		// Requirement: if observed is >= minObs && observed > expected then the
		// value is the difference, otherwise the value is zero.
		V observed = KmerHashTable<I, V>::values[y];

		if (observed >= minObs && observed > exp) {
			KmerHashTable<I, V>::values[y] = round(observed - exp);
		} else {
			KmerHashTable<I, V>::values[y] = 0;
		}
	}
}

} /* namespace nonltr */
