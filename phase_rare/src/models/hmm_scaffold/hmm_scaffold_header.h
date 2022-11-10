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

#ifndef _HMM_SCAFFOLD_H
#define _HMM_SCAFFOLD_H

#include <utils/otools.h>
#include <containers/conditioning_set/conditioning_set_header.h>
#include <containers/genotype_set/genotype_set_header.h>
#include <objects/hmm_parameters.h>
#include <containers/state_set.h>

#include <immintrin.h>

inline
float horizontal_add (const __m256& a) {
	__m128 vlow = _mm256_castps256_ps128(a);
	__m128 vhigh = _mm256_extractf128_ps(a, 1); // high 128
	vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
	__m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
	__m128 sums = _mm_add_ps(vlow, shuf);
	shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
	sums = _mm_add_ss(sums, shuf);    // (no wasted instructions, and all of them are the 4B minimum)
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
	std::vector < aligned_vector32 < float > > alpha;
	aligned_vector32 < float > beta;

public:
	//CONSTRUCTOR/DESTRUCTOR
	hmm_scaffold(variant_map & _V, genotype_set & _G, conditioning_set & _C, hmm_parameters & _M);
	~hmm_scaffold();

	void setup(unsigned int _hap);
	double forward();
	void backward(std::vector < std::vector < unsigned int > > & cevents, std::vector < int > & vpath);
	void viterbi(std::vector < int > & path);

};

#endif








