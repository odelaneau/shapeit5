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
#ifndef _GIBBS_SAMPLER_H
#define _GIBBS_SAMPLER__H

#include <utils/otools.h>
#include <containers/genotype_set.h>
#include <containers/state_set.h>
#include <containers/conditioning_set/conditioning_set_header.h>


class gibbs_sampler {
public:
	unsigned int nsamples;
	unsigned int niterations;
	unsigned int nburnin;
	bool rare;

	vector < bool > alleles;
	vector < bool > missing;
	//vector < bool > truth;

	vector < unsigned int > unphased;
	vector < vector < unsigned int > > cstates;
	vector < vector < float > > cprobs;
	vector < float > pprobs;

	gibbs_sampler(unsigned int, unsigned int, unsigned int);
	~gibbs_sampler();

	void loadRare(genotype_set & G, conditioning_set & C, state_set & P, unsigned int vr, float weight);
	void loadCommon(genotype_set & G, conditioning_set & C, state_set & P, unsigned int vc, float weight);

	void iterate(int &, int &);

	void pushRare(genotype_set & G, unsigned int vr);
	void pushCommon(genotype_set &, unsigned int vc);

};

#endif
