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
#include <objects/conditioning_set.h>

#define MAX_OVERLAP_HETS	0.9f
#define MIN_NUMBER_HETS		5

conditioning_set::conditioning_set(variant_map & _V, genotype_set & _G, haplotype_set & _H, unsigned int n_max_transitions, unsigned int n_max_missing) : V(_V), G(_G), H(_H) {
	phasing_probs = vector < double > (n_max_transitions, 0.0);
	missing_probs = vector < float > (n_max_missing * HAP_NUMBER, 0.0);

	O = vector < unsigned int > (H.n_total_haplotypes);
	iota(O.begin(), O.end(), 0);
	Oiter = 0;

	Svar = vector < vector < unsigned int > > (H.n_total_variants);
	Fvar = vector < bool > (H.n_total_variants);

}

conditioning_set::~conditioning_set() {
	free();
}

void conditioning_set::free () {
	vector < double > ().swap(phasing_probs);
	vector < float > ().swap(missing_probs);
	vector < unsigned int > ().swap(K);
	vector < unsigned int > ().swap(O);
}

void conditioning_set::make(unsigned int ind, bool process_rare, bool flip_rare) {
	//1. Update conditional haps
/*
	K.clear();
	unsigned long addr_offset = H.pbwt_nstored * H.n_target_samples * 2UL;
	vector < int > phap = vector < int > (2 * H.pbwt_depth, -1);
	for (int l = 0 ; l < H.pbwt_evaluated.size() ; l ++) {
		int abs_idx = H.pbwt_evaluated[l], rel_idx = H.pbwt_stored[l];
		if (rel_idx >= 0) {
			unsigned long curr_hap0 = 2*ind+0, curr_hap1 = 2*ind+1;
			for (int s = 0 ; s < H.pbwt_depth ; s ++) {
				int cond_hap0 = H.pbwt_neighbours[s * addr_offset + curr_hap0*H.pbwt_nstored + rel_idx];
				int cond_hap1 = H.pbwt_neighbours[s * addr_offset + curr_hap1*H.pbwt_nstored + rel_idx];
				if (cond_hap0 != phap[2*s+0]) { K.push_back(cond_hap0); phap[2*s+0] = cond_hap0; };
				if (cond_hap1 != phap[2*s+1]) { K.push_back(cond_hap1); phap[2*s+1] = cond_hap1; };
			}
		}
	}
*/
	K = H.pbwt_neighbours[ind];

	//2. Protect for IBD2
	ind_ibd2.clear();
	sort(K.begin(), K.end());
	K.erase(unique(K.begin(), K.end()), K.end());
	int count_het, match_het, nToBeRemoved = 0;
	vector < bool > vToBeRemoved = vector < bool > (K.size(), false);
	for (int k=1; k< K.size() ; k++) {
		unsigned int ind0 = K[k-1]/2;
		unsigned int ind1 = K[k]/2;
		if (ind0==ind1) {
			H.Hhap.getMatchHetCount(ind, ind0, count_het, match_het);
			assert(match_het <= count_het);
			float perc_matching_hets = (count_het - match_het) * 1.0f / count_het;
			if (perc_matching_hets > MAX_OVERLAP_HETS && count_het > MIN_NUMBER_HETS) {
				nToBeRemoved+=2;
				vToBeRemoved[k-1]=true;
				vToBeRemoved[k]=true;
				ind_ibd2.push_back(ind0);
			}
		}
		if (nToBeRemoved>0) {
			vector < unsigned int > Ktmp;
			Ktmp.reserve(K.size()-nToBeRemoved);
			for (int k=0; k< K.size() ; k++) if (!vToBeRemoved[k]) Ktmp.push_back(K[k]);
			K = Ktmp;
		}
	}

	//3. Add new random haps
	if (nToBeRemoved > 0) {
		for (int k = 0 ; k < nToBeRemoved ; k ++) {
			if ((O[Oiter]/2)!=ind) K.push_back(O[Oiter]);
			Oiter=(Oiter<(H.n_total_haplotypes-1))?Oiter+1:0;
		}
		sort(K.begin(), K.end());
		K.erase(unique(K.begin(), K.end()), K.end());
	}

	//4. Prepare Hhap and Hvar
	Hhap.reallocateFast(K.size(), H.n_common_variants);
	Hvar.reallocateFast(H.n_common_variants, K.size());
	Hhap.subset(H.Hhap, K);
	Hhap.transpose(Hvar);

	//5. Prepare Svar
	if (process_rare) {
		std::fill(Fvar.begin(), Fvar.end(), false);
		for (unsigned int vt = 0 ; vt < H.flag_common.size() ; vt ++) Svar[vt].clear();
		for (unsigned int vr = 0 ; vr < G.vecG[ind]->RareIndexes.size() ; vr ++) Fvar[G.vecG[ind]->RareIndexes[vr]] = true;
		if (!flip_rare) {
			for (unsigned int k = 0 ; k < K.size() ; k++) {
				for (unsigned int r = 0 ; r < H.Shap[K[k]].size() ; r++) {
					if (Fvar[H.Shap[K[k]][r]]) {
						Svar[H.Shap[K[k]][r]].push_back(k);
					}
				}
			}
		} else {
			for (unsigned int k = 0 ; k < K.size() ; k++) {
				//Pass 1 for main conditioning hap [usual shit]
				for (unsigned int r = 0 ; r < H.Shap[K[k]].size() ; r++) {
					if (Fvar[H.Shap[K[k]][r]]) {
						Svar[H.Shap[K[k]][r]].push_back(k);
					}
				}
				//Pass 2 for complementary conditioning hap [new shit] / Beagle trick that assumes IBD sharing for rare
				unsigned int complementary_k = MOD2(K[k])?(K[k]-1):(K[k]+1);
				for (unsigned int r = 0 ; r < H.Shap[complementary_k].size() ; r++) {
					if (Fvar[H.Shap[complementary_k][r]]) {
						Svar[H.Shap[complementary_k][r]].push_back(k);
					}
				}
			}
		}
	}
}

void conditioning_set::maskingTransitions(unsigned int ind, double error_rate) {
	vector < double > curr_transitions = vector < double > (4096, 0.0);
	unsigned int prev_dipcount = 1, curr_dipcount = 0, curr_transcount = 0;
	for (unsigned int s = 0, t = 0 ; s < G.vecG[ind]->n_segments ; s ++) {
		curr_dipcount = G.vecG[ind]->countDiplotypes(G.vecG[ind]->Diplotypes[s]);
		curr_transcount = prev_dipcount * curr_dipcount;

		double sumT = 0.0;
		for (unsigned int trel = 0 ; trel < curr_transcount ; trel ++) {
			curr_transitions[trel] = phasing_probs[t+trel] * (G.vecG[ind]->ProbabilityMask[t+trel]?(1.0-error_rate):(error_rate));
			sumT += curr_transitions[trel];
		}

		if (sumT > numeric_limits<double>::min() && !isnan(sumT))
			for (unsigned int trel = 0 ; trel < curr_transcount ; trel ++)
				phasing_probs[t+trel] = curr_transitions[trel] / sumT;

		t += curr_transcount;
		prev_dipcount = curr_dipcount;
	}
}
