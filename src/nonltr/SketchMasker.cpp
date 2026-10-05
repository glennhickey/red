/*
 * SketchMasker.cpp
 *
 * See SketchMasker.h.  The counting and chaining follow the skc2 prototype exactly, so for the same seeds,
 * thresholds and table size the masked intervals are the same.
 */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>

#include "SketchMasker.h"
#include "ChromListMaker.h"
#include "../utility/Location.h"
#include "../exception/InvalidStateException.h"

using namespace std;
using namespace utility;
using namespace exception;

namespace nonltr {

namespace {

const int BATCH = 8192;      // table lookups batched for prefetching
const int PREFETCH = 24;     // how far ahead in a batch to prefetch
const long BLOCK = 1L << 22; // positions scored at a time
const double REFERENCE_SIZE = 5.923e9; // genome size the thresholds were tuned at (sawshark)

// A/C/G/T in either case -> 0..3, anything else -> 4 (breaks the seed)
struct Encoding {
	uint8_t code[256];
	Encoding() {
		memset(code, 4, sizeof code);
		code[(int) 'A'] = code[(int) 'a'] = 0;
		code[(int) 'C'] = code[(int) 'c'] = 1;
		code[(int) 'G'] = code[(int) 'g'] = 2;
		code[(int) 'T'] = code[(int) 't'] = 3;
	}
};
const Encoding ENC;

inline uint64_t mix64(uint64_t x) {
	x ^= x >> 30;
	x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27;
	x *= 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

// The low bit of every 2-bit base, packed together: A=0 C=1 G=0 T=1, i.e. purine 0, pyrimidine 1.
// Equivalent to BMI2's pext(x, 0x5555555555555555) without needing it.
inline uint64_t evenBits(uint64_t x) {
	x &= 0x5555555555555555ULL;
	x = (x | (x >> 1)) & 0x3333333333333333ULL;
	x = (x | (x >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
	x = (x | (x >> 4)) & 0x00FF00FF00FF00FFULL;
	x = (x | (x >> 8)) & 0x0000FFFF0000FFFFULL;
	x = (x | (x >> 16)) & 0x00000000FFFFFFFFULL;
	return x;
}

// Near-periodic seeds (period 1 to 4, in bases, and for the pattern also in purine/pyrimidine space) are
// skipped when counting and scoring: they are low-complexity sequence, which Red handles on its own.
// Thresholds are the fraction of positions that differ: about 75% for random bases, 50% for random R/Y.
inline bool lowComplexity(uint64_t fwd, int span, bool ry) {
	for (int per = 1; per <= 4; per++) {
		int n = span - per;
		uint64_t m = (n >= 32) ? ~0ULL : ((1ULL << (2 * n)) - 1);
		uint64_t d = (fwd ^ (fwd >> (2 * per))) & m;
		int nd = __builtin_popcountll((d | (d >> 1)) & 0x5555555555555555ULL & m);
		if (nd * 100 < 35 * n) {
			return true;
		}
	}
	if (ry) {
		uint64_t r = evenBits(fwd);
		for (int per = 1; per <= 4; per++) {
			int n = span - per;
			uint64_t m = (1ULL << n) - 1;
			int nd = __builtin_popcountll((r ^ (r >> per)) & m);
			if (nd * 100 < 22 * n) {
				return true;
			}
		}
	}
	return false;
}

// Union of two sorted, disjoint interval lists
void mergeInto(vector<pair<int, int> > & a, const vector<pair<int, int> > & b) {
	if (b.empty()) {
		return;
	}
	vector<pair<int, int> > out;
	out.reserve(a.size() + b.size());
	size_t i = 0, j = 0;
	while (i < a.size() || j < b.size()) {
		const pair<int, int> & next = (j == b.size() || (i < a.size() && a[i].first <= b[j].first)) ? a[i++] : b[j++];
		if (!out.empty() && next.first <= out.back().second) {
			out.back().second = max(out.back().second, next.second);
		} else {
			out.push_back(next);
		}
	}
	a.swap(out);
}

}

SketchMasker::SketchMasker(const vector<string> & fileList) {
	cout << endl << endl;
	cout << "Stage 0: Sketch ..." << endl;

	// Genome size, and the number of records in each file
	genomeSize = 0;
	string header;
	string seq;
	bool hadSequence = false;
	regions.resize(fileList.size());
	for (size_t f = 0; f < fileList.size(); f++) {
		ChromListMaker maker(fileList.at(f));
		while (maker.nextSequence(header, seq, hadSequence)) {
			genomeSize += (long) seq.size();
			regions[f].push_back(vector<pair<int, int> >());
		}
	}
	if (genomeSize == 0) {
		cout << "Sketch: no sequence, nothing to mask." << endl;
		return;
	}

	// The table grows with the genome up to 2 x 2^30 cells; past that, a fixed table is enough because
	// the thresholds grow with the genome too, so the collision noise stays the same fraction of them.
	int bits = (int) ceil(log2((double) genomeSize)) + 1;
	bits = max(20, min(30, bits));

	const Seed seeds[2] = {
		{ "16-mer", 16, false, 128.0, 16, 20, 3 },
		{ "32-base purine/pyrimidine pattern", 32, true, 64.0, 8, 20, 6 } };
	for (int i = 0; i < 2; i++) {
		runSeed(seeds[i], fileList, bits);
	}
}

void SketchMasker::runSeed(const Seed & seed, const vector<string> & fileList, int bits) {
	const int K = seed.span;
	long c = lround(seed.scale * (double) genomeSize / REFERENCE_SIZE);
	const uint32_t C = (uint32_t) min(65535L, max((long) seed.floor, c));
	const int width = C <= 255 ? 1 : 2;
	const uint64_t M = 1ULL << bits;
	const uint64_t MM = M - 1;
	const uint64_t spanMask = K == 32 ? ~0ULL : ((1ULL << (2 * K)) - 1);

	cout << "Sketch: " << seed.name << ", genome " << genomeSize << " bp, threshold " << C
			<< ", table 2x2^" << bits << " cells of " << width << " byte(s)" << endl;

	void * table = calloc((size_t) M * 2, width);
	if (table == NULL) {
		throw InvalidStateException(string("Could not allocate the sketch table."));
	}
	uint8_t * a8 = (uint8_t *) table;
	uint8_t * b8 = a8 + M;
	uint16_t * a16 = (uint16_t *) table;
	uint16_t * b16 = a16 + M;

	vector<uint32_t> i1(BATCH), i2(BATCH), ip(BATCH);
	string header;
	string seq;
	bool hadSequence = false;

	// ---- count every seed of the genome ----
	for (size_t f = 0; f < fileList.size(); f++) {
		ChromListMaker maker(fileList.at(f));
		while (maker.nextSequence(header, seq, hadSequence)) {
			const char * s = seq.data();
			const long n = (long) seq.size();
			uint64_t fwd = 0, rc = 0;
			int run = 0, nb = 0;
			for (long j = 0; j <= n; j++) {
				if (j < n) {
					uint8_t b = ENC.code[(unsigned char) s[j]];
					if (b > 3) {
						run = 0;
						fwd = rc = 0;
						continue;
					}
					fwd = ((fwd << 2) | b) & spanMask;
					rc = (rc >> 2) | ((uint64_t) (b ^ 3) << (2 * (K - 1)));
					if (++run < K || lowComplexity(fwd, K, seed.ry)) {
						continue;
					}
					uint64_t kf = seed.ry ? evenBits(fwd) : fwd;
					uint64_t kr = seed.ry ? evenBits(rc) : rc;
					uint64_t h = mix64(kf < kr ? kf : kr);
					i1[nb] = (uint32_t) (h & MM);
					i2[nb] = (uint32_t) ((h >> 32) & MM);
					nb++;
				}
				if (nb == BATCH || (j == n && nb > 0)) {
					if (width == 1) {
						for (int q = 0; q < nb; q++) {
							if (q + PREFETCH < nb) {
								__builtin_prefetch(a8 + i1[q + PREFETCH], 1);
								__builtin_prefetch(b8 + i2[q + PREFETCH], 1);
							}
							if (a8[i1[q]] < 255) a8[i1[q]]++;
							if (b8[i2[q]] < 255) b8[i2[q]]++;
						}
					} else {
						for (int q = 0; q < nb; q++) {
							if (q + PREFETCH < nb) {
								__builtin_prefetch(a16 + i1[q + PREFETCH], 1);
								__builtin_prefetch(b16 + i2[q + PREFETCH], 1);
							}
							if (a16[i1[q]] < 65535) a16[i1[q]]++;
							if (b16[i2[q]] < 65535) b16[i2[q]]++;
						}
					}
					nb = 0;
				}
			}
		}
	}

	// ---- score every record: chain the seeds whose count reaches the threshold ----
	vector<uint32_t> cb(BLOCK);
	long masked = 0, added = 0;
	for (size_t f = 0; f < fileList.size(); f++) {
		ChromListMaker maker(fileList.at(f));
		size_t r = 0;
		while (maker.nextSequence(header, seq, hadSequence)) {
			const char * s = seq.data();
			const long n = (long) seq.size();
			vector<pair<int, int> > found;
			long cs = -1, ce = -1;
			int cn = 0;
			uint64_t fwd = 0, rc = 0;
			int run = 0, nb = 0;
			for (long b0 = 0; b0 < n; b0 += BLOCK) {
				const long len = min(BLOCK, n - b0);
				for (long j = 0; j < len; j++) {
					cb[j] = 0;
					uint8_t b = ENC.code[(unsigned char) s[b0 + j]];
					if (b > 3) {
						run = 0;
						fwd = rc = 0;
						continue;
					}
					fwd = ((fwd << 2) | b) & spanMask;
					rc = (rc >> 2) | ((uint64_t) (b ^ 3) << (2 * (K - 1)));
					if (++run < K || lowComplexity(fwd, K, seed.ry)) {
						continue;
					}
					uint64_t kf = seed.ry ? evenBits(fwd) : fwd;
					uint64_t kr = seed.ry ? evenBits(rc) : rc;
					uint64_t h = mix64(kf < kr ? kf : kr);
					i1[nb] = (uint32_t) (h & MM);
					i2[nb] = (uint32_t) ((h >> 32) & MM);
					ip[nb] = (uint32_t) j;
					nb++;
					if (nb == BATCH || j == len - 1) {
						for (int q = 0; q < nb; q++) {
							uint32_t x, y;
							if (width == 1) {
								if (q + PREFETCH < nb) {
									__builtin_prefetch(a8 + i1[q + PREFETCH], 0);
									__builtin_prefetch(b8 + i2[q + PREFETCH], 0);
								}
								x = a8[i1[q]];
								y = b8[i2[q]];
							} else {
								if (q + PREFETCH < nb) {
									__builtin_prefetch(a16 + i1[q + PREFETCH], 0);
									__builtin_prefetch(b16 + i2[q + PREFETCH], 0);
								}
								x = a16[i1[q]];
								y = b16[i2[q]];
							}
							cb[ip[q]] = x < y ? x : y;
						}
						nb = 0;
					}
				}
				// seeds whose last base came before a trailing run of non-ACGT
				for (int q = 0; q < nb; q++) {
					uint32_t x = width == 1 ? a8[i1[q]] : a16[i1[q]];
					uint32_t y = width == 1 ? b8[i2[q]] : b16[i2[q]];
					cb[ip[q]] = x < y ? x : y;
				}
				nb = 0;
				// cb[j] is the count of the seed starting at b0 + j - (K - 1)
				for (long j = 0; j < len; j++) {
					const long S = b0 + j - (K - 1);
					if (S < 0) {
						continue;
					}
					const bool hi = cb[j] >= C;
					if (hi && cs >= 0 && S <= ce + seed.gap) {
						if (S + K > ce) {
							ce = S + K;
						}
						cn++;
						continue;
					}
					if (cs >= 0 && hi) {
						if (cn >= seed.minHigh) {
							found.push_back(make_pair((int) cs, (int) min(ce, n)));
						}
						cs = -1;
					}
					if (hi) {
						cs = S;
						ce = S + K;
						cn = 1;
					}
				}
			}
			if (cs >= 0 && cn >= seed.minHigh) {
				found.push_back(make_pair((int) cs, (int) min(ce, n)));
			}

			vector<pair<int, int> > & kept = regions[f][r];
			long before = 0, after = 0;
			for (size_t q = 0; q < kept.size(); q++) {
				before += kept[q].second - kept[q].first;
			}
			for (size_t q = 0; q < found.size(); q++) {
				masked += found[q].second - found[q].first;
			}
			mergeInto(kept, found);
			for (size_t q = 0; q < kept.size(); q++) {
				after += kept[q].second - kept[q].first;
			}
			added += after - before;
			r++;
		}
	}
	free(table);

	cout << "Sketch: " << seed.name << " masks " << masked << " bp (" << (100.0 * masked / genomeSize)
			<< "% of the genome), " << added << " bp of it not masked by an earlier seed" << endl;
}

SketchMasker::~SketchMasker() {
}

void SketchMasker::getRegions(int f, int r, vector<ILocation *> * out) const {
	if (f < 0 || f >= (int) regions.size() || r < 0 || r >= (int) regions[f].size()) {
		return;
	}
	const vector<pair<int, int> > & v = regions[f][r];
	for (size_t q = 0; q < v.size(); q++) {
		out->push_back(new Location((int) v[q].first, (int) (v[q].second - 1)));
	}
}

long SketchMasker::getGenomeSize() const {
	return genomeSize;
}

} /* namespace nonltr */
