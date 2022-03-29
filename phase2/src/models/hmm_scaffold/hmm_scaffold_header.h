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
#ifndef _HMM_SCAFFOLD_H
#define _HMM_SCAFFOLD_H

#include <utils/otools.h>
#include <containers/conditioning_set/conditioning_set_header.h>
#include <containers/genotype_set.h>
#include <objects/hmm_parameters.h>
#include <containers/state_set.h>
#include <boost/align/aligned_allocator.hpp>
#include <immintrin.h>


template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

inline
float horizontal_add (const __m256& a) {
	__m128 vlow = _mm256_castps256_ps128(a);
	__m128 vhigh = _mm256_extractf128_ps(a, 1); // high 128
	vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
	__m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
	__m128 sums = _mm_add_ps(vlow, shuf);
	shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
	sums = _mm_add_ss(sums, shuf);    // (no wasted instructions, and all of them are the 4B minimum)

	/*
	aligned_vector32 < float > tmp  = aligned_vector32 < float > (8);
	_mm256_store_ps(&tmp[0], a);
	cout << _mm_cvtss_f32(sums) << " " << tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7] << endl;
	 */

	return _mm_cvtss_f32(sums);
}

class hmm_scaffold {
public:
	//DATA
	conditioning_set & C;
	genotype_set & G;
	variant_map & V;
	hmm_parameters & M;

	//CONSTANT
	unsigned int hap;
	float match_prob[2];
	unsigned int nstates;

	//
	bitmatrix Hvar;
	bitmatrix Hhap;

	//ARRAYS
	vector < aligned_vector32 < float > > alpha;
	aligned_vector32 < float > beta;
	vector < bool > storageEvents;


public:
	//CONSTRUCTOR/DESTRUCTOR
	hmm_scaffold(unsigned int _hap, variant_map & _V, genotype_set & _G, conditioning_set & _C, hmm_parameters & _M);
	~hmm_scaffold();

	void forward();
	void backward(vector < bool > & cevents, vector < state > & cstates);
};

#endif








