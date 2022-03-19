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

#include <models/gibbs_sampler/gibbs_sampler_header.h>

class hmm_scaffold {
public:
	//EXTERNAL DATA
	conditioning_set & C;
	genotype_set & G;
	variant_map & V;
	hmm_parameters & M;

	//GIBBS SAMPLER
	gibbs_sampler MCMC;

	//CONSTANT
	float emit[2];
	float rev_emit[2];

	//ARRAYS
	vector < bool > bufferA0;
	vector < bool > bufferA1;
	vector < unsigned int > K;
	vector < unsigned long > nstates;
	vector < unsigned long > sstates;
	vector < float > alpha;
	vector < float > beta;

public:
	//CONSTRUCTOR/DESTRUCTOR
	hmm_scaffold(variant_map & _V, genotype_set & _G, conditioning_set & _C, hmm_parameters & _M);
	~hmm_scaffold();

	//
	void prefetchAlleles0(unsigned int vs);
	void prefetchAlleles1(unsigned int vs);

	void forward();
	void backward(float threshold, unsigned int niterations, unsigned int nburnin);

	//
	void forward_initTransitions(unsigned int h);
	void forward_updateTransitions(unsigned int vs, unsigned int h);
	void forward_updateEmission(unsigned int vs, unsigned int h);
	void forward_normalize(unsigned int vs, unsigned int h);
	void forward_reverseEmission(unsigned int vs, unsigned int h);
	void forward_reverseTransitions(unsigned int vs, unsigned int h);
	void backward_initTransitions(unsigned int h);
	void backward_updateTransitions(unsigned int vs, unsigned int h);
	void backward_updateEmission(unsigned int vs, unsigned int h);
	void backward_normalize(unsigned int vs, unsigned int h);
	void getAlphaBetaProduct(unsigned int h, vector < float > & alphaXbeta);

};

inline
void hmm_scaffold::prefetchAlleles0(unsigned int vs) {
	for (int h = 0 ; h < C.n_haplotypes ; h ++) bufferA0[h] = C.Hvar.get(vs, h);
}

inline
void hmm_scaffold::prefetchAlleles1(unsigned int vs) {
	for (int h = 0 ; h < C.n_haplotypes ; h ++) bufferA1[h] = C.Hvar.get(vs, h);
}

#endif








