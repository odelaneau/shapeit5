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

hmm_scaffold::hmm_scaffold(unsigned int _hap, variant_map & _V, genotype_set & _G, conditioning_set & _C, hmm_parameters & _M) : V(_V), G(_G), C(_C), M(_M){
	hap = _hap;
	emit[0] = 1.0f; emit[1] = M.ed/M.ee;

	unsigned max_nstates = 0;
	for (int h = 0 ; h < C.n_haplotypes ; h ++) if (C.indexes_pbwt_neighbour[h].size() > max_nstates) max_nstates = C.indexes_pbwt_neighbour[h].size();

	alpha = vector < vector < float > > (C.n_scaffold_variants, vector < float > (max_nstates, 0.0f));
	alphaSum = vector < float > (C.n_scaffold_variants, 0.0f);
	beta = vector < float > (max_nstates, 1.0f);

	nstates = C.indexes_pbwt_neighbour[hap].size();
	Hhap.subset(C.Hhap, C.indexes_pbwt_neighbour[hap]);
	Hvar.allocate(C.n_scaffold_variants, nstates);
	Hhap.transpose(Hvar);
}

hmm_scaffold::~hmm_scaffold() {
	alpha.clear();
	alphaSum.clear();
	beta.clear();
	storageEvents.clear();
}

void hmm_scaffold::forward() {
	for (int vs = 0 ; vs < C.n_scaffold_variants ; vs ++) {
		if (!vs) fill(alpha[vs].begin(), alpha[vs].begin() + nstates, 1.0f / nstates);
		else {
			float f0 = M.t[vs-1] / nstates;
			float  f1 = M.nt[vs-1] / alphaSum[vs-1];
			for (int k = 0 ; k < nstates ; k ++) alpha[vs][k] = alpha[vs-1][k] * f1 + f0;
		}
		alphaSum[vs] = 0.0f;
		for (int k = 0 ; k < nstates ; k ++) {
			alpha[vs][k] *= emit[C.Hhap.get(hap, vs) != Hvar.get(vs, k)];
			alphaSum[vs] += alpha[vs][k];
		}
	}
}

unsigned int hmm_scaffold::backward(vector < bool > & cevents, vector < cprobs > & cstates) {
	unsigned int n_total_states = 0;
	float sum = 0.0f, scale = 0.0f, threshold = 1.0f / nstates;
	vector < float > alphaXbeta_curr = vector < float >(nstates, 0.0f);
	vector < float > alphaXbeta_prev = vector < float >(nstates, 0.0f);

	cstates.clear();
	fill (beta.begin(), beta.end(), 1.0f / nstates);

	for (int vs = C.n_scaffold_variants - 1 ; vs >= 0 ; vs --) {

		//Transitions
		if (vs < C.n_scaffold_variants - 1) {
			float f0 = M.t[vs] / nstates;
			float  f1 = M.nt[vs] / sum;
			for (int k = 0 ; k < nstates ; k ++) beta[k] = beta[k] * f1 + f0;
		}

		//Products
		scale = 0.0f;
		for (int k = 0 ; k < nstates ; k ++) {
			alphaXbeta_curr[k] = alpha[vs][k] * beta[k];
			scale += alphaXbeta_curr[k];
		}
		scale = 1.0f / scale;
		for (int k = 0 ; k < nstates ; k ++) {
			alphaXbeta_curr[k] *= scale;
			//cout << hap << " " << vs << " " << k << " " << alphaXbeta_curr[k] << endl;
		}

		//Emission
		sum = 0.0f;
		for (int k = 0 ; k < nstates ; k ++) {
			beta[k] *= emit[C.Hhap.get(hap, vs) != Hvar.get(vs, k)];
			sum += beta[k];
		}

		//Storage
		if (cevents[vs+1]) {
			if (vs == C.n_scaffold_variants-1) copy(alphaXbeta_curr.begin(), alphaXbeta_curr.begin() + nstates, alphaXbeta_prev.begin());
			unsigned int n_reserved_states = 0;
			for (int k = 0 ; k < nstates ; k ++) if (alphaXbeta_curr[k] >= threshold || alphaXbeta_prev[k] >= threshold) n_reserved_states ++;
			cstates.emplace_back(n_reserved_states, hap, vs+2);
			for (int k = 0 ; k < nstates ; k ++)
				if (alphaXbeta_curr[k] >= threshold || alphaXbeta_prev[k] >= threshold)
					cstates.back().push(k, (unsigned char)(alphaXbeta_curr[k] * 254), (unsigned char)(alphaXbeta_prev[k] * 254));
			n_total_states += cstates.back().size();
		}

		//Saving products
		copy(alphaXbeta_curr.begin(), alphaXbeta_curr.begin() + nstates, alphaXbeta_prev.begin());
	}

	if (cevents[0]) {
		unsigned int n_reserved_states = 0;
		for (int k = 0 ; k < nstates ; k ++) if (alphaXbeta_curr[k] >= threshold) n_reserved_states ++;
		cstates.emplace_back(n_reserved_states, hap, 1);
		for (int k = 0 ; k < nstates ; k ++)
			if (alphaXbeta_curr[k] >= threshold)
				cstates.back().push(k, (unsigned char)(alphaXbeta_curr[k] * 254), (unsigned char)(alphaXbeta_curr[k] * 254));
		n_total_states += cstates.back().size();
	}

	return n_total_states;
}
