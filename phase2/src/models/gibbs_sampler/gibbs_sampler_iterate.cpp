////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#include <models/gibbs_sampler/gibbs_sampler_header.h>

void gibbs_sampler::iterate(float weight) {
	float prob_hap0 [2], prob_hap1 [2], sum;
	vector < float > prob_genotype = vector < float > (4);


	for (int iter = 0 ; iter < niterations ; iter ++) {

		for (int ui = 0 ; ui < unphased_indexes.size() ; ui ++) {
			unsigned int ind = unphased_indexes[ui];

			//Sanity check debugging
			assert(missing [ind] || (alleles0[ind] != alleles1[ind]));

			//Compute imputation probabilities
			prob_hap0[0] = prob_hap0[1] = 0.0f;
			for (int k = 0 ; k < state_indexes[2*ind+0].size() ; k ++) {
				prob_hap0[alleles0[state_indexes[2*ind+0][k]]] += state_lprobs[2*ind+0][k] * weight + state_rprobs[2*ind+0][k] * (1.0f - weight);
				prob_hap1[alleles1[state_indexes[2*ind+1][k]]] += state_lprobs[2*ind+1][k] * weight + state_rprobs[2*ind+1][k] * (1.0f - weight);
			}

			//Samples phased genotype
			fill(prob_genotype.begin(), prob_genotype.end(), 1.0f);
			if (!missing[ind]) {
				prob_genotype[0] = 0.0f;
				prob_genotype[3] = 0.0f;
			}

			prob_genotype[0] *= prob_hap0[0] * prob_hap1[0];
			prob_genotype[1] *= prob_hap0[0] * prob_hap1[1];
			prob_genotype[2] *= prob_hap0[1] * prob_hap1[0];
			prob_genotype[3] *= prob_hap0[1] * prob_hap1[1];
			sum = accumulate(prob_genotype.begin(), prob_genotype.end(), 0.0f);
			assert(sum > 0);

			switch (rng.sample(prob_genotype, sum)) {
			case 0:	alleles0[ind] = 0; alleles1[ind] = 0; break;
			case 1:	alleles0[ind] = 0; alleles1[ind] = 1; break;
			case 2:	alleles0[ind] = 1; alleles1[ind] = 0; break;
			case 3:	alleles0[ind] = 1; alleles1[ind] = 1; break;
			}

			if (iter >= nburnin) {
				phasing_probs[4*ind+0] += prob_genotype[0];
				phasing_probs[4*ind+1] += prob_genotype[1];
				phasing_probs[4*ind+2] += prob_genotype[2];
				phasing_probs[4*ind+3] += prob_genotype[3];
			}
		}
	}

	for (int ui = 0 ; ui < unphased_indexes.size() ; ui ++) {
		unsigned int ind = unphased_indexes[ui];
		unsigned int bestg = distance(phasing_probs.begin() + 4*ind, max_element(phasing_probs.begin() + 4*ind, phasing_probs.begin() + 4*(ind+1))) % 4;
		switch (bestg) {
		case 0:	alleles0[ind] = 0; alleles1[ind] = 0; break;
		case 1:	alleles0[ind] = 0; alleles1[ind] = 1; break;
		case 2:	alleles0[ind] = 1; alleles1[ind] = 0; break;
		case 3:	alleles0[ind] = 1; alleles1[ind] = 1; break;
		}
	}
}
