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
#include <models/haplotype_segment.h>


haplotype_segment::haplotype_segment(genotype * _G, conditioning_set & _C, hmm_parameters & _M) : G(_G), C(_C), M(_M){
	n_states = C.K.size();
	n_missing = G->n_missing;
	n_rarevar = G->RareIndexes.size();
	n_segment = G->n_segments;

	probSumT = 0.0f;
#ifdef __AVX2__
	prob = aligned_vector32 < float > (HAP_NUMBER * n_states, 0.0f);
	probSumH = aligned_vector32 < float > (HAP_NUMBER, 0.0f);
	probSumK = aligned_vector32 < float > (n_states, 0.0f);
	Alpha = vector < aligned_vector32 < float > > (n_segment, aligned_vector32 < float > (HAP_NUMBER * n_states, 0.0f));
	AlphaSum = vector < aligned_vector32 < float > > (n_segment, aligned_vector32 < float > (HAP_NUMBER, 0.0f));
	AlphaSumSum = aligned_vector32 < float > (n_segment, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER * n_states, 0.0f));
		AlphaSumMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER, 0.0f));
	}
#else
	prob = vector < float > (HAP_NUMBER * n_states, 0.0f);
	probSumH = vector < float > (HAP_NUMBER, 0.0f);
	probSumK = vector < float > (n_states, 0.0f);
	Alpha = vector < vector < float > > (n_segment, vector < float > (HAP_NUMBER * n_states, 0.0f));
	AlphaSum = vector < vector < float > > (n_segment, vector < float > (HAP_NUMBER, 0.0f));
	AlphaSumSum = vector < float > (n_segment, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < vector < float > > (n_missing, vector < float > (HAP_NUMBER * n_states, 0.0f));
		AlphaSumMissing = vector < vector < float > > (n_missing, vector < float > (HAP_NUMBER, 0.0f));
	}
#endif
	if (n_rarevar > 0) {
		AlphaRare = vector < vector < float > > (n_rarevar, vector < float > (2 * n_states, 1.0f));
	}
}

haplotype_segment::~haplotype_segment() {
	G = NULL;

	n_states = 0;
	n_segment = 0;
	n_missing = 0;
	n_rarevar = 0;

	prob.clear();
	probSumK.clear();
	probSumH.clear();
	Alpha.clear();
	AlphaSum.clear();
	AlphaSumSum.clear();
	AlphaMissing.clear();
	AlphaSumMissing.clear();
	AlphaRare.clear();
}

void haplotype_segment::forwardRephase() {
	float yt, nt, f1, f20, f21, cm, mm = M.ed/M.ee;
	float emit[2]; emit[0] = 1.0f; emit[1] = mm;

	for (int vt = 0, vc = 0, vr = 0 ; vt < C.H.flag_common.size() ; vt ++) {

		if (C.H.flag_common[vt]) {
			bool ag0 = VAR_GET_HAP0(MOD2(vc), G->Variants[DIV2(vc)]);
			bool ag1 = VAR_GET_HAP1(MOD2(vc), G->Variants[DIV2(vc)]);

			if (vc == 0) {
				fill (prob.begin(), prob.begin() + 2*n_states, 1.0f / n_states);
				probSum0 = probSum1 = 1.0f;
			} else {
				f1 = M.t[vc-1] / n_states;
				f20 = M.nt[vc-1] / probSum0;
				f21 = M.nt[vc-1] / probSum1;
				probSum0 = probSum1 = 0.0f;
				for (int k = 0 ; k < n_states ; k ++) {
					bool ah = C.Hvar.get(vc, k);
					prob[2*k+0] = (prob[2*k+0] * f20 + f1) * emit[ah!=ag0];
					prob[2*k+1] = (prob[2*k+1] * f21 + f1) * emit[ah!=ag1];
					probSum0 += prob[2*k+0];
					probSum1 += prob[2*k+1];
				}
			}
			vc ++;
		} else if (G->RareIndexes[vr] == vt) {
			if (C.Svar[vt].size() > 0) {

				// ----- IMPORTANT -----
				sort(C.Svar[vt].begin(), C.Svar[vt].end());
				C.Svar[vt].erase(unique(C.Svar[vt].begin(), C.Svar[vt].end()), C.Svar[vt].end());
				// ----- ---------------

				bool het = VAR_GET_HET(MOD2(vr), G->RareVariants[DIV2(vr)]);
				bool mis = VAR_GET_MIS(MOD2(vr), G->RareVariants[DIV2(vr)]);
				if (het || mis) {
					if (vc != 0) {
						//cout << vt << " " << vc - 1 << " " << C.V.vec_full[vt]->cm << " " << C.V.vec_common[vc-1]->cm << endl;
						cm = C.V.vec_full[vt]->cm - C.V.vec_common[vc-1]->cm;
						if (cm == 0.0f) cm = 1e-7;
						assert(cm > 0.0f);
						yt = M.getTransProb(cm);
						nt = 1.0f - yt;
						f1 = yt / n_states;
						f20 = nt / probSum0;
						f21 = nt / probSum1;
						for (int k = 0 ; k < n_states ; k ++) {
							AlphaRare[vr][2*k+0] *= (prob[2*k+0] * f20 + f1);
							AlphaRare[vr][2*k+1] *= (prob[2*k+1] * f21 + f1);
						}
					}
				}
			}
			vr ++;
		}
	}
}

void haplotype_segment::backwardRephase(bool maximize) {
	float yt, nt, f1, f20, f21, cm, sc, mm = M.ed/M.ee;
	float emit[2]; emit[0] = 1.0f; emit[1] = mm;
	double probs [2][2];

	for (int vt = C.H.flag_common.size() - 1, vc = G->n_common_variants - 1, vr = G->RareIndexes.size() - 1 ; vt >= 0 ; vt --) {

		if (C.H.flag_common[vt]) {
			bool ag0 = VAR_GET_HAP0(MOD2(vc), G->Variants[DIV2(vc)]);
			bool ag1 = VAR_GET_HAP1(MOD2(vc), G->Variants[DIV2(vc)]);

			if (vc == (G->n_common_variants - 1)) {
				fill (prob.begin(), prob.begin() + 2*n_states, 1.0f / n_states);
				probSum0 = probSum1 = 1.0f;
			} else {		//AVX2 this thing
				f1 = M.t[vc] / n_states;
				f20 = M.nt[vc] / probSum0;
				f21 = M.nt[vc] / probSum1;
				probSum0 = probSum1 = 0.0f;
				for (int k = 0 ; k < n_states ; k ++) {
					bool ah = C.Hvar.get(vc, k);
					prob[2*k+0] = (prob[2*k+0] * f20 + f1) * emit[ah!=ag0];
					prob[2*k+1] = (prob[2*k+1] * f21 + f1) * emit[ah!=ag1];
					probSum0 += prob[2*k+0];
					probSum1 += prob[2*k+1];
				}
			}
			vc --;
		} else if (G->RareIndexes[vr] == vt) {
			bool het = VAR_GET_HET(MOD2(vr), G->RareVariants[DIV2(vr)]);
			bool mis = VAR_GET_MIS(MOD2(vr), G->RareVariants[DIV2(vr)]);
			bool minor = C.V.vec_full[vt]->minor;
			if (het || mis) {
				probs[0][0] = probs[0][1] = probs[1][0] = probs[1][1] = 0.5f;

				if (C.Svar[vt].size() > 0) {
					vector < float > probs2 = vector < float > (2 * n_states, 1.0f);
					float psum0 = 1.0f, psum1 = 1.0f;

					if (vc < (G->n_common_variants - 1)) {
						psum0 = psum1 = 0.0f;
						cm = C.V.vec_common[vc+1]->cm - C.V.vec_full[vt]->cm;
						if (cm == 0.0f) cm = 1e-7;
						assert(cm > 0.0f);
						yt = M.getTransProb(cm);
						nt = 1.0f - yt;
						f1 = yt / n_states;
						f20 = nt / probSum0;
						f21 = nt / probSum1;
						for (int k = 0 ; k < n_states ; k ++) {
							probs2[2*k+0] *= (prob[2*k+0] * f20 + f1);
							probs2[2*k+1] *= (prob[2*k+1] * f21 + f1);
							psum0 += probs2[2*k+0];
							psum1 += probs2[2*k+1];
						}
					}

					psum0 = 1.0f / psum0;
					psum1 = 1.0f / psum1;
					assert(!isnan(psum0));
					assert(!isnan(psum1));

					for (int k = 0 ; k < n_states ; k ++) {
						//cout << k << " " << AlphaRare[vr][2*k+0] << " " << AlphaRare[vr][2*k+1] << " " << probs2[2*k+0] << " " << probs2[2*k+1] <<endl;
						AlphaRare[vr][2*k+0] *= (probs2[2*k+0] * psum0);
						AlphaRare[vr][2*k+1] *= (probs2[2*k+1] * psum1);
						//cout << k << " " << AlphaRare[vr][2*k+0] << " " << AlphaRare[vr][2*k+1] << " " << probs2[2*k+0] << " " << probs2[2*k+1] <<endl;
					}
					probs[0][0] = probs[0][1] = probs[1][0] = probs[1][1] = 0.0f;
					//This loop assumes C.Svar[vt] to be sorted!
					for (int k = 0, r = 0 ; k < n_states ; k ++) {
						bool which_allele = ((r < C.Svar[vt].size()) && (C.Svar[vt][r] == k)) ^ (!minor);
						probs[0][which_allele] += AlphaRare[vr][2*k+0];
						probs[1][which_allele] += AlphaRare[vr][2*k+1];
						r += ((r < C.Svar[vt].size()) && (C.Svar[vt][r] == k));
					}
				}
				vector < double > genotype_probabilities = vector < double > (4, 1.0f);
				if (het) {	//Ensuring R|R or A|A not sampled for hets
					genotype_probabilities[0] = 0.0f;
					genotype_probabilities[3] = 0.0f;
				}
				genotype_probabilities[0] *= probs[0][0] * probs[1][0];
				genotype_probabilities[1] *= probs[0][0] * probs[1][1];
				genotype_probabilities[2] *= probs[0][1] * probs[1][0];
				genotype_probabilities[3] *= probs[0][1] * probs[1][1];
				double sum_genotype_probabilities =  genotype_probabilities[0] + genotype_probabilities[1] + genotype_probabilities[2] + genotype_probabilities[3];

				/*
				if (C.V.vec_full[vt]->getMAC() > 1) {
					cout << het << " " << mis << " " << C.V.vec_full[vt]->getMAC()
						<< " " << stb.str(genotype_probabilities[0]/sum_genotype_probabilities, 3)
						<< " " << stb.str(genotype_probabilities[1]/sum_genotype_probabilities, 3)
						<< " " << stb.str(genotype_probabilities[2]/sum_genotype_probabilities, 3)
						<< " " << stb.str(genotype_probabilities[3]/sum_genotype_probabilities, 3) << endl;
				}
				assert(!isnan(sum_genotype_probabilities));
				assert(sum_genotype_probabilities>0.0f);
				assert(het && ((genotype_probabilities[1] + genotype_probabilities[2])>=0));
				*/

				unsigned int genotype;
				if (!maximize) genotype = rng.sample(genotype_probabilities, sum_genotype_probabilities);
				else genotype = std::max_element(genotype_probabilities.begin(), genotype_probabilities.end()) - genotype_probabilities.begin();

				switch (genotype) {
				case 0:	VAR_CLR_HAP0(MOD2(vr), G->RareVariants[DIV2(vr)]);
						VAR_CLR_HAP1(MOD2(vr), G->RareVariants[DIV2(vr)]);
						break;
				case 1:	VAR_CLR_HAP0(MOD2(vr), G->RareVariants[DIV2(vr)]);
						VAR_SET_HAP1(MOD2(vr), G->RareVariants[DIV2(vr)]);
						break;
				case 2:	VAR_SET_HAP0(MOD2(vr), G->RareVariants[DIV2(vr)]);
						VAR_CLR_HAP1(MOD2(vr), G->RareVariants[DIV2(vr)]);
						break;
				case 3:	VAR_SET_HAP0(MOD2(vr), G->RareVariants[DIV2(vr)]);
						VAR_SET_HAP1(MOD2(vr), G->RareVariants[DIV2(vr)]);
						break;
				}
			}
			vr --;
		}
	}
}


void haplotype_segment::forward() {
	forward_pass = true;
	curr_segment_index = 0;
	curr_segment_locus = 0;
	curr_ambiguous = 0;
	curr_missing = 0;

	for (curr_locus = 0 ; curr_locus < G->n_common_variants ; curr_locus++) {
		bool amb = VAR_GET_AMB(MOD2(curr_locus), G->Variants[DIV2(curr_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_locus), G->Variants[DIV2(curr_locus)]);
		bool hom = !(amb || mis);

		if (curr_locus == 0) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != 0) {
			if (hom) RUN_HOM();
			else if (amb) RUN_AMB();
			else RUN_MIS();
		} else {
			if (hom) COLLAPSE_HOM();
			else if (amb) COLLAPSE_AMB();
			else  COLLAPSE_MIS();
		}

		if (curr_segment_locus == (G->Lengths[curr_segment_index] - 1)) SUMK();
		if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) {
			Alpha[curr_segment_index] = prob;
			AlphaSum[curr_segment_index] = probSumH;
			AlphaSumSum[curr_segment_index] = probSumT;
		}
		if (mis) {
			AlphaMissing[curr_missing] = prob;
			AlphaSumMissing[curr_missing] = probSumH;
			curr_missing ++;
		}

		curr_segment_locus ++;
		curr_ambiguous += amb;
		if (curr_segment_locus >= G->Lengths[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
}

int haplotype_segment::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities) {
	forward_pass = false;
	int n_underflow_recovered = 0;
	curr_segment_index = n_segment - 1;
	curr_segment_locus = G->Lengths.back() - 1;
	curr_ambiguous = G->n_ambiguous - 1;
	curr_missing = G->n_missing - 1;
	curr_transition = G->n_transitions - 1;

	for (curr_locus = G->n_common_variants - 1 ; curr_locus >= 0 ; curr_locus--) {
		bool amb = VAR_GET_AMB(MOD2(curr_locus), G->Variants[DIV2(curr_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_locus), G->Variants[DIV2(curr_locus)]);
		bool hom = !(amb || mis);

		if (curr_locus == (G->n_common_variants - 1)) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != G->Lengths[curr_segment_index] - 1) {
			if (hom) RUN_HOM();
			else if (amb) RUN_AMB();
			else RUN_MIS();
		} else {
			if (hom) COLLAPSE_HOM();
			else if (amb) COLLAPSE_AMB();
			else COLLAPSE_MIS();
		}
		if (curr_segment_locus == 0) SUMK();
		if (curr_locus == 0) SET_FIRST_TRANS(transition_probabilities);

		if (curr_segment_locus == 0 && curr_locus != 0) {
			int ret = SET_OTHER_TRANS(transition_probabilities);
			if (ret < 0) return ret;
			else n_underflow_recovered += ret;
		}

		if (mis) {
			IMPUTE(missing_probabilities);
			curr_missing--;
		}


		curr_segment_locus--;
		curr_ambiguous -= amb;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths[curr_segment_index] - 1;
		}
	}
	return n_underflow_recovered;
}

void haplotype_segment::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
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

int haplotype_segment::SET_OTHER_TRANS(vector < double > & transition_probabilities) {
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
	curr_transition -= (n_transitions - 1);
	for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_transition + t] = (DProbs[t] * scaleDip);
	curr_transition --;
	return underflow_recovered;
}
