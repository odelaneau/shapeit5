/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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

#include <containers/haplotype_set.h>
#include <containers/genotype_set.h>
#include <containers/variant_map.h>

class conditioning_set {
public:
	//Static data structures
	variant_map & V;
	genotype_set & G;
	haplotype_set & H;

	//Dynamic data structures
	bitmatrix Hhap, Hvar;
	vector < vector < unsigned int > > Svar;
	vector < bool > Fvar;

	//Probability vectors
	vector < double > phasing_probs;
	vector < float > missing_probs;

	//Conditioning states
	vector < unsigned int > K;
	vector < unsigned int > O;
	int Oiter;

	//IBD2 tracking
	vector < int > ind_ibd2;

	conditioning_set(variant_map & , genotype_set & , haplotype_set & , unsigned int n_max_transitions , unsigned int n_max_missing);
	~conditioning_set();

	void free();

	void make(unsigned int, bool, bool);
	void maskingTransitions(unsigned int, double);
};

#endif
