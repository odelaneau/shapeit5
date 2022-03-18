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
#include <models/hmm_scaffold/hmm_scaffold_header.h>

hmm_scaffold::hmm_scaffold(variant_map & _V, genotype_set * _G, conditioning_set & _C, hmm_parameters & _M) : V(_V), G(_G), C(_C), M(_M){
	//
	nstates = vector < unsigned long > (C.n_haplotypes, 0);
	sstates = vector < unsigned long > (C.n_haplotypes, 0);
	for (int i = 0; i < G.n_samples ; i ++) nstates[i] = C.indexes_pbwt_neighbour[i].size();
	for (int i = 1; i < G.n_samples ; i ++) sstates[i] = sstates[i-1] + C.indexes_pbwt_neighbour[i-1].size();
	n_total_states = sstates.back() + nstates.back();

	//
	emit[0] = 1.0f;
	emit[1] = M.ed/M.ee;
	rev_emit[0] = 1.0f;
	rev_emit[1] = M.ee/M.ed;

	//
	K = vector < unsigned int > (n_total_states);
	for (unsigned int h = 0, k = 0 ; h < C.n_haplotypes ; h ++) {
		copy(K.begin() + k, C.indexes_pbwt_neighbour[h/2].begin(), C.indexes_pbwt_neighbour[h/2].end());
		k += C.indexes_pbwt_neighbour[h/2].size();
	}

	//
	bufferA0 = vector < bool > (C.n_haplotypes, false);
	bufferA1 = vector < bool > (C.n_haplotypes, false);
	alpha = vector < float > (n_total_states, 0.0f);
	beta = vector < float > (n_total_states, 1.0f);
	cprobs = vector < vector < prob_state > > (C.n_haplotypes);
}

hmm_scaffold::~hmm_scaffold() {
	nstates.clear();
	sstates.clear();
	alpha.clear();
	beta.clear();
}

void hmm_scaffold::forward() {
	tac.clock();
	vrb.title("Forward pass");

	//Main loop
	for (int vs = 0 ; vs < C.n_scaffold_variants ; vs ++) {

		//Prefetch alleles at scaffold to speed up compute
		prefetchAlleles0(vs);

		//
		for (int h = 0 ; h < C.n_haplotypes ; h ++) {

			//Update forward with transition from l-1 to l
			if (vs) forward_updateTransitions(vs, h);
			else forward_initTransitions(h);

			//Update forward with emission at l
			forward_updateEmission(vs, h);

			//Normalize forward
			forward_normalize(vs, h);
		}
		vrb.progress("  * Processing", (vs+1)*1.0/C.n_scaffold_variants);
	}
	vrb.bullet("Timing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void hmm_scaffold::backward(float threshold) {
	tac.clock();

	//Variables
	vector < unsigned int > VC, VR;
	vector < float > alphaXbeta_prev, alphaXbeta_curr;

	vrb.title("Backward pass");

	//Initialization of betas
	for (int h = 0 ; h < C.n_haplotypes ; h ++) backward_initTransitions(h);

	//Main loop
	for (int vs = C.n_scaffold_variants - 2 ; vs >= 0 ; vs --) {

		//Prefetch alleles at scaffold to speed up compute
		prefetchAlleles0(vs);
		prefetchAlleles1(vs+1);

		//Get rare and common unphased variants in the [l, l+1] interval
		V.getCommonVariants(vs, vs+1, VC);
		V.getRareVariants(vs, vs+1, VR);

		//Are there any M/m or ./. at these variants
		G.getUnphasedIndexes(VC, VR, IDX);

		//If yes prepare MCMC data structures


		//
		for (int h = 0, i = 0 ; h < C.n_haplotypes ; h ++) {

			//If there are unphased genotypes in [l, l+1] for sample h/2
			if (IDX.size() > 0 && IDX[i] == h/2) {
				//Allocate vector to store alphaXbeta products
				alphaXbeta_prev = vector < float > (nstates[h], 0.0f);
				alphaXbeta_curr = vector < float > (nstates[h], 0.0f);

				//Get alphaXbeta product at l+1 [redundant / to be optimized out]
				getAlphaBetaProduct(alphaXbeta_prev);
			}

			//Update backward with emission at l+1
			backward_updateEmission(vs+1, h);

			//Update backward with transition from l+1 to l
			backward_updateTransitions(vs, h);

			//Normalize backward
			backward_normalize(vs, h);

			//Reverse forward from l+1 to l
			forward_reverseEmission(vs, h);
			forward_reverseTransitions(vs, h);
			forward_normalize(vs, h);

			//If there are unphased genotypes in [l, l+1] for sample h/2
			if (IDX.size() > 0 && IDX[i] == h/2) {
				//Get alphaXbeta product at l
				getAlphaBetaProduct(alphaXbeta_curr);

				//Compress alphaXbeta products at l and l + 1

			}
		}

		vrb.progress("  * Processing", (C.n_scaffold_variants-vs)*1.0/C.n_scaffold_variants);
	}
	vrb.bullet("Timing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}








			//Compression
			for (int k = 0 ; k < nstates[h] ; k ++)
				if (alphaXbeta_prev[k] >= threshold || alphaXbeta_curr[k] >= threshold)
					cprobs[h].emplace_back(k, alphaXbeta_curr[k], alphaXbeta_prev[k]);











		}
	}
}










			}





		}


		bool a1 = C.HSvar.get(vs, 2 * i_sample + 1);
		if (vs == 0) {
			fill (alpha.begin(), alpha.begin() + 2 * n_states, 1.0f / n_states);
		} else {
			for (int k = 0 ; k < n_states ; k ++) {
				bool ah = C.Hvar.get(vc, k);
				alpha[(2*vs*n_states) + 2*k + 0] = (alpha[(2*(vs-1)*n_states) + 2*k + 0] * f20 + f1) * emit[ah!=a0];
				alpha[(2*vs*n_states) + 2*k + 1] = (alpha[(2*(vs-1)*n_states) + 2*k + 1] * f20 + f1) * emit[ah!=a1];
				alphaSum[2*(vs-1) + 0] += alpha[(2*vs*n_states) + 2*k + 0];
				alphaSum[2*(vs-1) + 1] += alpha[(2*vs*n_states) + 2*k + 1];
			}
		}
	}
}

void hmm_scaffold::backward(vector < cfloat > & compressed_probabilities_hap0, vector < cfloat > & compressed_probabilities_hap1, float threshold) {
	float emit [2], betaSum_prev[2], betaSum_curr[2], prodSum [2], scale [2];
	emit[0] = 1.0f;
	emit[1] = M.ed/M.ee;
	vector < float > beta = vector < float > (2 * n_states, 0.0f);
	vector < float > alphaXbeta_curr = vector < float > (2 * n_states, 0.0f);
	vector < float > alphaXbeta_prev = vector < float > (2 * n_states, 0.0f);
	for (int vs = C.n_scaffold_variants - 1 ; vs >= 0 ; vs --) {
		bool a0 = C.HSvar.get(vs, 2 * i_sample + 0);
		bool a1 = C.HSvar.get(vs, 2 * i_sample + 1);

		//Backward recursion part 1
		if (vs == (C.n_scaffold_variants - 1)) {
			fill (beta.begin(), beta.begin() + 2 * n_states, 1.0f / n_states);
			betaSum_curr[0] = betaSum_curr[1] = 1.0f;
		} else {
			float f1 = M.t[vs] / n_states;
			float f20 = M.nt[vs] / betaSum_prev[0];
			float f21 = M.nt[vs] / betaSum_prev[1];
			for (int k = 0 ; k < n_states ; k ++) {
				beta[2*k + 0] = ([2*k + 0] * f20 + f1);
				beta[2*k + 1] = ([2*k + 1] * f20 + f1);
			}
		}

		//Product of forward and backward probs
		prodSum[0] = prodSum[1] = 1.0f;
		scale[0] = 1.0f * alphaSum[2*vs+0];
		scale[1] = 1.0f * alphaSum[2*vs+1];
		for (int k = 0 ; k < 2*n_states ; k ++) {
			alphaXbeta_curr[k] = (alpha[(2*vs*n_states) + k] * scale[k%2]) * beta[k];
			prodSum[k%2] += alphaXbeta_curr[k];
		}

		//Normalization of product
		prodSum[0] = 1.0f / prodSum[0];
		prodSum[1] = 1.0f / prodSum[1];
		for (int k = 0 ; k < 2*n_states ; k ++) alphaXbeta_curr[k] *= prodSum[k%2];

		//Compression
		if (vs < (C.n_scaffold_variants - 1)) {
			for (int k = 0 ; k < n_states ; k ++) {
				if ((alphaXbeta_curr[2*k+0] >= threshold) || (alphaXbeta_prev[2*k+0] >= threshold)) {

				}
			}
		}

		//Backward recursion part 2
		betaSum_curr[0] = betaSum_curr[1] = 0.0f;
		for (int k = 0 ; k < n_states ; k ++) {
			bool ah = C.Hvar.get(vc, k);
			beta[2*k + 0] *= emit[ah!=a0];
			beta[2*k + 1] *= emit[ah!=a1];
			betaSum_curr[0] += beta[2*k + 0];
			betaSum_curr[1] += beta[2*k + 1];
		}

		//Saving beta sums
		betaSum_prev[0] = betaSum_curr[0];
		betaSum_prev[1] = betaSum_curr[1];
		alphaXbeta_prev = alphaXbeta_curr;
	}
}
