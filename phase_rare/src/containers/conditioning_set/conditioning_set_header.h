/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _CONDITIONING_SET_H
#define _CONDITIONING_SET_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>
#include <containers/genotype_set/genotype_set_header.h>

class cflip {
public:
	unsigned int pgenotype;
	float support;

	cflip() {
		pgenotype = 0;
		support = 0.0f;
	}

	cflip(unsigned int _pgenotype, float _support) {
		pgenotype = _pgenotype;
		support = _support;
	}

	void set(unsigned int _pgenotype, float _support) {
		pgenotype = _pgenotype;
		support = _support;
	}

	bool betterThan(cflip & cf) {
		if (abs(support) > abs(cf.support)) return true;
		else return false;
	}
};


class conditioning_set : public haplotype_set {
public:
	//VARIANT INDEXING
	std::vector < bool > sites_pbwt_evaluation;
	std::vector < bool > sites_pbwt_selection;
	std::vector < int > sites_pbwt_grouping;
	unsigned int sites_pbwt_ngroups;
	unsigned long int ncollisions;
	unsigned long int npushes;

	//PARAMETERS FOR PBWT
	int depth_common;
	int depth_rare;

	//PHASE DATA
	std::vector < cflip > CF;

	//STATE DATA
	unsigned int shuffledI;
	std::vector < unsigned int > shuffledO;
	std::vector < std::pair < unsigned int, unsigned int > > indexes_pbwt_neighbour_serialized;
	std::vector < std::vector < unsigned int > > indexes_pbwt_neighbour;

	//IBD2 protection
	std::vector < std::vector < unsigned int > > IBD2;

	//CONSTRUCTOR/DESTRUCTOR
	conditioning_set();
	~conditioning_set();
	void initialize(variant_map &, float, float, int, int, int);

	//STATES PROCESSING
	void storeCommon(std::vector < int > & A, std::vector < int > & M);
	void storeRare(std::vector < int > & R, std::vector < sparse_genotype > & G);
	void select(variant_map &, genotype_set & G);

	//IBD2 scanning
	void scanIBD2(variant_map &);
	bool checkIBD2(unsigned int, unsigned int);

	//SOLVING
	/*
	void solveRare1(std::vector < int > &, std::vector < int > &, genotype_set &, unsigned int);
	void solveRare2(std::vector < int > &, std::vector < int > &, std::vector < int > &, genotype_set &, unsigned int, unsigned int, std::vector < float > &);
	*/
	void solveRareForward(std::vector < int > & A, std::vector < int > & C, std::vector < int > & R, genotype_set & G, unsigned int vr, float vr_cm, std::vector < float > & vs_cm);
	void solveRareBackward(std::vector < int > & A, std::vector < int > & C, std::vector < int > & R, genotype_set & G, unsigned int vr, float vr_cm, std::vector < float > & vs_cm);
	void solve(variant_map &, genotype_set &);
};

inline
bool conditioning_set::checkIBD2(unsigned int h1, unsigned int h2) {
	int i1 = std::min(h1/2, h2/2);
	int i2 = std::max(h1/2, h2/2);
	if (i1 == i2) return true;
	for (int i  = 0 ; i < IBD2[i1].size() ; i ++) if (IBD2[i1][i] == i2) return true;
	return false;
}

#endif
