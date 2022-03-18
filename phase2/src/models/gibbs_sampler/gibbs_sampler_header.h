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
#define _GIBBS_SAMPLER_H

#include <containers/genotype_set.h>

class gibbs_sampler {
public:
	//MCMC parameters
	unsigned int nsamples;
	unsigned int niterations;
	unsigned int nburnin;

	//Sample indexes to be phased
	vector < unsigned int > unphased_indexes;			// #unphased

	//Workspace for computation
	vector < bool > alleles0;							// #haplotypes / Phased alleles for unphased genotypes
	vector < bool > alleles1;							// #haplotypes / Phased alleles for unphased genotypes
	vector < bool > missing;

	//Compressed state space w/ copying probs
	vector < vector < unsigned int > > state_indexes;	// #samples / conditioning states / indexes
	vector < vector < float > > state_lprobs;			// #samples / conditioning states / left probs
	vector < vector < float > > state_rprobs;			// #samples / conditioning states / right probs

	//Phasing probs
	vector < float > phasing_probs;

	//CONSTRUCTOR/DESTRUCTOR
	gibbs_sampler(unsigned int);
	~gibbs_sampler();

	//INPUT
	bool loadCommonUnphasedGenotypes(unsigned int, genotype_set &);
	bool loadRareUnphasedGenotypes(unsigned int, genotype_set &, bool);
	void loadStateSpace(vector < vector < unsigned int > > &, vector < vector < float > > &, vector < vector < float > > &, float);

	//MCMC
	void iterate();

	//OUTPUT
	void pushCommonPhasedGenotypes(unsigned int, genotype_set &);
	void pushRarePhasedGenotypes(unsigned int, genotype_set &);
};

#endif
