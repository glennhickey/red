/*
 * ChromosomeOneDigit.cpp
 *
 *  Created on: Jul 31, 2012
 *      Author: Hani Zakaria Girgis, PhD at the NCB1/NLM/NIH
 * A	A
 * T	T
 * G	G
 * C	C
 * R	G or A
 * Y	T or C
 * M	A or C
 * K	G or T
 * S	G or C
 * W	A or T
 * H	A or C or T
 * B	G or T or C
 * V	G or C or A
 * D	G or T or A
 * N	G or T or A or C
 */
#include <iostream>
#include <map>

#include "Chromosome.h"
#include "ChromosomeOneDigit.h"
#include "../exception/InvalidInputException.h"

using namespace exception;

namespace nonltr {

ChromosomeOneDigit::ChromosomeOneDigit() :
		Chromosome() {
}

ChromosomeOneDigit::ChromosomeOneDigit(string fileName) :
		Chromosome(fileName) {
	help();
}

ChromosomeOneDigit::ChromosomeOneDigit(string seq, string info) :
		Chromosome(seq, info) {
	help();
}

void ChromosomeOneDigit::help() {
	// Build codes
	buildCodes();
	// Modify the sequence in the super class
	encodeNucleotides();
}

void ChromosomeOneDigit::finalize() {
	Chromosome::finalize();
	help();
}

/**
 * A 256-entry view of the "codes" map, indexed by the byte itself.  The map
 * lookups it replaces cost two red-black-tree walks per base (count() then
 * at()), which on a multi-gigabase genome is the single most expensive thing
 * in the parse.  INVALID_CODE marks a byte with no entry in the map, so the
 * same InvalidInputException is still thrown for it.
 */
const char ChromosomeOneDigit::INVALID_CODE = (char) -1;

void ChromosomeOneDigit::buildCodeTable() {
	for (int i = 0; i < 256; i++) {
		codeTable[i] = INVALID_CODE;
	}
	for (map<char, char>::const_iterator it = codes->begin();
			it != codes->end(); it++) {
		codeTable[(unsigned char) it->first] = it->second;
	}
}

void ChromosomeOneDigit::buildCodes() {
	// Make map
	codes = new map<char, char>();

	// Certain nucleotides
	codes->insert(map<char, char>::value_type('A', (char) 0));
	codes->insert(map<char, char>::value_type('C', (char) 1));
	codes->insert(map<char, char>::value_type('G', (char) 2));
	codes->insert(map<char, char>::value_type('T', (char) 3));

	// Common uncertain nucleotide
	// codes->insert(map<char, char>::value_type('N', (char) 4));

	// Uncertain nucleotides
	codes->insert(map<char, char>::value_type('R', codes->at('G')));
	codes->insert(map<char, char>::value_type('Y', codes->at('C')));
	codes->insert(map<char, char>::value_type('M', codes->at('A')));
	codes->insert(map<char, char>::value_type('K', codes->at('T')));
	codes->insert(map<char, char>::value_type('S', codes->at('G')));
	codes->insert(map<char, char>::value_type('W', codes->at('T')));
	codes->insert(map<char, char>::value_type('H', codes->at('C')));
	codes->insert(map<char, char>::value_type('B', codes->at('T')));
	codes->insert(map<char, char>::value_type('V', codes->at('A')));
	codes->insert(map<char, char>::value_type('D', codes->at('T')));
	codes->insert(map<char, char>::value_type('N', codes->at('C')));
	codes->insert(map<char, char>::value_type('X', codes->at('G')));

	buildCodeTable();
}

ChromosomeOneDigit::~ChromosomeOneDigit() {
	codes->clear();
	delete codes;
}

/**
 * This method converts nucleotides in the segments to single digit codes
 */
void ChromosomeOneDigit::encodeNucleotides() {

  char * b = &base[0];

  for (int s = 0; s < segment->size(); s++) {
    int segStart = segment->at(s)->at(0);
    int segEnd = segment->at(s)->at(1);
    for (int i = segStart; i <= segEnd; i++) {
      char code = codeTable[(unsigned char) b[i]];
      if (code != INVALID_CODE) {
	b[i] = code;
      } else {
	string msg = "Invalid nucleotide: ";
	msg.append(1, b[i]);
	throw InvalidInputException(msg);
      }
    }
  }

  // Digitize skipped segments
  int segNum = segment->size();
  if(segNum > 0){
    // The first interval - before the first segment
    int segStart = 0; 
    int segEnd = segment->at(0)->at(0)-1; 

    for (int s = 0; s <= segNum; s++) {      
      for (int i = segStart; i <= segEnd; i++) {
	char c = b[i];
	if(c != 'N'){
	  char code = codeTable[(unsigned char) c];
	  if (code != INVALID_CODE) {
	    b[i] = code;
	  } else {
	    string msg = "Invalid nucleotide: ";
	    msg.append(1, c);
	    throw InvalidInputException(msg);
	  }
	}
      }

      // The regular intervals between two segments
      if(s < segNum-1){
	segStart = segment->at(s)->at(1)+1;
	segEnd = segment->at(s+1)->at(0)-1;
      }
      // The last interval - after the last segment
      else if(s == segNum - 1){
	segStart = segment->at(s)->at(1)+1;
	segEnd = base.size()-1;
      } 
    } 
  }
}

/*
void ChromosomeOneDigit::encodeNucleotides() {
	int seqLen = base.size();

	for (int i = 0; i < seqLen; i++) {
		if (codes->count(base[i]) > 0) {
			base[i] = codes->at(base[i]);
		} else {
			string msg = "Invalid nucleotide: ";
			msg.append(1, base[i]);
			throw InvalidInputException(msg);
		}
	}

}
*/

/**
 * Cannot be called on already finalized object.
 */
void ChromosomeOneDigit::makeR() {
	//cout << "Making reverse ..." << endl;
	makeReverse();
	reverseSegments();
}

/**
 * Cannot be called on already finalized object.
 */
void ChromosomeOneDigit::makeRC() {
	//cout << "Making reverse complement ..." << endl;
	makeComplement();
	makeReverse();
	reverseSegments();
}

void ChromosomeOneDigit::makeComplement() {
	// A byte table rather than a map, for the same reason as codeTable: this
	// runs once per base per strand, twice per sequence in the scan stage.
	char complement[256];
	for (int i = 0; i < 256; i++) {
		complement[i] = INVALID_CODE;
	}

	// Certain nucleotides
	complement[0] = (char) 3;
	complement[1] = (char) 2;
	complement[2] = (char) 1;
	complement[3] = (char) 0;

	// Unknown nucleotide
	complement[(unsigned char) 'N'] = 'N';

	// Convert a sequence to its complement
	size_t seqLen = base.size();
	char * b = &base[0];
	for (size_t i = 0; i < seqLen; i++) {
		char c = complement[(unsigned char) b[i]];
		if (c != INVALID_CODE) {
			b[i] = c;
		} else {
			cerr << "Error: The digit " << (char) b[i];
			cerr << " does not represent a base." << endl;
			exit(2);
		}
	}
}

void ChromosomeOneDigit::makeReverse() {
	size_t last = base.size() - 1;

	// Last index to be switched
	size_t middle = base.size() / 2;

	char * b = &base[0];
	for (size_t i = 0; i < middle; i++) {
		char temp = b[last - i];
		b[last - i] = b[i];
		b[i] = temp;
	}
}

void ChromosomeOneDigit::reverseSegments() {
	int segNum = segment->size();
	int lastBase = size() - 1;

	// Calculate the coordinate on the main strand
	for (int i = 0; i < segNum; i++) {
		vector<int> * seg = segment->at(i);

		int s = lastBase - seg->at(1);
		int e = lastBase - seg->at(0);
		seg->clear();
		seg->push_back(s);
		seg->push_back(e);
	}

	// Reverse the regions within the list
	int lastRegion = segNum - 1;
	int middle = segNum / 2;
	for (int i = 0; i < middle; i++) {
		vector<int> * temp = segment->at(lastRegion - i);
		(*segment)[lastRegion - i] = segment->at(i);
		(*segment)[i] = temp;
	}
}

}
