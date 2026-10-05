/*
 * SketchMasker.h
 *
 * -sketch: also mask sequence whose seeds are high-copy across the whole genome, judged from their counts
 * alone, so no alignment is needed.  Two seeds are counted into count-min sketches over every record:
 * a 16-mer, and a 32-base pattern of purines and pyrimidines, which still matches old repeat copies that
 * have drifted by transitions.  A position whose seed count reaches the threshold is "high"; high seeds
 * closer than 20 bp are chained, and chains with enough high seeds (3 for the 16-mer, 6 for the pattern)
 * are masked.
 *
 * The thresholds scale with genome size, which keeps them a fixed density per lastz chunk: for a genome of
 * G bases, max(16, round(128 G / 5.923e9)) for the 16-mer and max(8, round(64 G / 5.923e9)) for the
 * pattern.  The floors keep small genomes from masking every duplicated k-mer.  Cells are one byte, or two
 * once a threshold is above 255 (genomes over about 12 Gbp), so the threshold can always be reached.
 *
 * All of this runs before Red builds its own table, one seed at a time, so the sketch's table is freed
 * before Red's is allocated.  What is kept is the masked intervals of every record.
 */

#ifndef SKETCHMASKER_H_
#define SKETCHMASKER_H_

#include <string>
#include <vector>
#include <utility>
#include <cstdint>

#include "../utility/ILocation.h"

using namespace std;
using namespace utility;

namespace nonltr {

class SketchMasker {
private:
	// [file][record] -> half-open intervals, sorted and disjoint.  int, like Red's own coordinates.
	vector<vector<vector<pair<int, int> > > > regions;
	long genomeSize;

	struct Seed {
		const char * name;
		int span;       // bases covered by one seed
		bool ry;        // purine/pyrimidine pattern instead of the full bases
		double scale;   // threshold per 5.923 Gbp of genome
		int floor;      // lowest threshold
		int gap;        // high seeds closer than this are chained
		int minHigh;    // chains with fewer high seeds are not masked
	};

	void runSeed(const Seed &, const vector<string> &, int);

public:
	SketchMasker(const vector<string> & fileList);
	virtual ~SketchMasker();

	// Appends the masked regions of record r of file f to out, as Locations with inclusive ends, sorted
	// and disjoint.  Files and records are numbered in the order fileList and the files list them.
	void getRegions(int f, int r, vector<ILocation *> * out) const;

	long getGenomeSize() const;
};

} /* namespace nonltr */

#endif /* SKETCHMASKER_H_ */
