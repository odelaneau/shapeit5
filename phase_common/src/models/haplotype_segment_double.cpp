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

#include <models/haplotype_segment_double.h>

using namespace std;

haplotype_segment_double::haplotype_segment_double(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, window & W, hmm_parameters & _M) : G(_G), M(_M){
	segment_first = W.start_segment;
	segment_last = W.stop_segment;
	locus_first = W.start_locus;
	locus_last = W.stop_locus;
	ambiguous_first = W.start_ambiguous;
	ambiguous_last = W.stop_ambiguous;
	missing_first = W.start_missing;
	missing_last = W.stop_missing;
	transition_first = W.start_transition;
	transition_last = W.stop_transition;
	n_cond_haps = idxH.size();
	n_missing = missing_last - missing_first + 1;

	probSumT = 0.0f;
	prob = aligned_vector32 < double > (HAP_NUMBER * n_cond_haps, 0.0f);
	probSumH = aligned_vector32 < double > (HAP_NUMBER, 0.0f);
	probSumK = aligned_vector32 < double > (n_cond_haps, 0.0f);
	Alpha = vector < aligned_vector32 < double > > (segment_last - segment_first + 1, aligned_vector32 < double > (HAP_NUMBER * n_cond_haps, 0.0f));
	AlphaLocus = vector < int > (segment_last - segment_first + 1, 0);
	AlphaSum = vector < aligned_vector32 < double > > (segment_last - segment_first + 1, aligned_vector32 < double > (HAP_NUMBER, 0.0f));
	AlphaSumSum = aligned_vector32 < double > (segment_last - segment_first + 1, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < aligned_vector32 < double > > (n_missing, aligned_vector32 < double > (HAP_NUMBER * n_cond_haps, 0.0f));
		AlphaSumMissing = vector < aligned_vector32 < double > > (n_missing, aligned_vector32 < double > (HAP_NUMBER, 0.0f));
	}
	//Cache efficient data transfer for conditioning haplotypes
	curr_rel_locus_offset = Hhap.subset(H, idxH, locus_first, locus_last);
	Hvar.allocateFast(Hhap.n_cols, Hhap.n_rows);
	Hhap.transpose(Hvar);
}

haplotype_segment_double::~haplotype_segment_double() {
	G = NULL;
	segment_first = 0;
	segment_last = 0;
	locus_first = 0;
	locus_last = 0;
	ambiguous_first = 0;
	ambiguous_last = 0;
	transition_first = 0;
	n_cond_haps = 0;
	n_missing = 0;
	curr_segment_index = 0;
	curr_segment_locus = 0;
	curr_abs_locus = 0;
	curr_rel_locus = 0;
	curr_abs_ambiguous = 0;
	curr_abs_transition = 0;
	probSumT = 0.0;
	prob.clear();
	probSumK.clear();
	probSumH.clear();
	Alpha.clear();
	AlphaSum.clear();
	AlphaSumSum.clear();
	AlphaMissing.clear();
	AlphaSumMissing.clear();
}

void haplotype_segment_double::forward() {
	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	curr_abs_ambiguous = ambiguous_first;
	curr_abs_missing = missing_first;
	prev_abs_locus = locus_first;

	for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool update_prev_locus = true;
		char rare_allele = M.rare_allele[curr_abs_locus];
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool hom = !(amb || mis);
		yt = (curr_abs_locus == locus_first)?0.0:M.getForwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;

		if (curr_rel_locus == 0) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != 0) {
			if (hom) update_prev_locus = RUN_HOM(rare_allele);
			else if (amb) RUN_AMB();
			else RUN_MIS();
		} else {
			if (hom) COLLAPSE_HOM();
			else if (amb) COLLAPSE_AMB();
			else  COLLAPSE_MIS();
		}
		prev_abs_locus=update_prev_locus?curr_abs_locus:prev_abs_locus;

		if (curr_segment_locus == (G->Lengths[curr_segment_index] - 1)) SUMK();
		if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) {
			Alpha[curr_segment_index - segment_first] = prob;
			AlphaSum[curr_segment_index - segment_first] = probSumH;
			AlphaSumSum[curr_segment_index - segment_first] = probSumT;
			AlphaLocus[curr_segment_index - segment_first] = prev_abs_locus;
		}
		if (mis) {
			AlphaMissing[curr_rel_missing] = prob;
			AlphaSumMissing[curr_rel_missing] = probSumH;
			curr_abs_missing ++;
		}

		curr_segment_locus ++;
		curr_abs_ambiguous += amb;
		if (curr_segment_locus >= G->Lengths[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
}

int haplotype_segment_double::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities) {
	int n_underflow_recovered = 0;
	curr_segment_index = segment_last;
	curr_segment_locus = G->Lengths[segment_last] - 1;
	curr_abs_ambiguous = ambiguous_last;
	curr_abs_missing = missing_last;
	curr_abs_transition = transition_last;
	prev_abs_locus = locus_last;

	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		char rare_allele = M.rare_allele[curr_abs_locus];
		bool update_prev_locus = true;
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool hom = !(amb || mis);
		yt = (curr_abs_locus == locus_last)?0.0:M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;

		if (curr_abs_locus == locus_last) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != G->Lengths[curr_segment_index] - 1) {
			if (hom) update_prev_locus = RUN_HOM(rare_allele);
			else if (amb) RUN_AMB();
			else RUN_MIS();
		} else {
			if (hom) COLLAPSE_HOM();
			else if (amb) COLLAPSE_AMB();
			else COLLAPSE_MIS();
		}
		if (curr_segment_locus == 0) SUMK();
		prev_abs_locus=update_prev_locus?curr_abs_locus:prev_abs_locus;

		if (curr_abs_locus == 0) SET_FIRST_TRANS(transition_probabilities);
		if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
			int ret = SET_OTHER_TRANS(transition_probabilities);
			if (ret < 0) return ret;
			else n_underflow_recovered += ret;
		}

		if (mis) {
			IMPUTE(missing_probabilities);
			curr_abs_missing--;
		}


		curr_segment_locus--;
		curr_abs_ambiguous -= amb;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths[curr_segment_index] - 1;
		}
	}
	return n_underflow_recovered;
}

void haplotype_segment_double::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	double scale = 1.0f / probSumT, scaleDip = 0.0f;
	unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	vector < double > cprobs = vector < double > (n_transitions, 0.0);
	for (unsigned int d = 0, t = 0 ; d < 64 ; ++d) {
		if (DIP_GET(G->Diplotypes[0], d)) {
			cprobs[t] = (double)(probSumH[DIP_HAP0(d)]*scale) * (double)(probSumH[DIP_HAP1(d)]*scale);
			scaleDip += cprobs[t++];
		}
	}
	scaleDip = 1.0f / scaleDip;
	for (unsigned int t = 0 ; t < n_transitions ; t ++) transition_probabilities[t] = cprobs[t] * scaleDip;
}

int haplotype_segment_double::SET_OTHER_TRANS(vector < double > & transition_probabilities) {
	int underflow_recovered = 0;
	if (TRANS_HAP()) return -1;
	if (TRANS_DIP_MULT()) {
		if (TRANS_DIP_ADD()) return -2;
		else underflow_recovered = 1;
	}
	unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
	unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
	unsigned int n_transitions = curr_dipcount * prev_dipcount;
	double scaleDip = 1.0 / sumDProbs;
	curr_abs_transition -= (n_transitions - 1);
	for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_abs_transition + t] = DProbs[t] * scaleDip;
	curr_abs_transition --;
	return underflow_recovered;
}
