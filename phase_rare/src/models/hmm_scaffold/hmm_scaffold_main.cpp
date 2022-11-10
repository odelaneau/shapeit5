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

#include <models/hmm_scaffold/hmm_scaffold_header.h>

using namespace std;

hmm_scaffold::hmm_scaffold(variant_map & _V, genotype_set & _G, conditioning_set & _C, hmm_parameters & _M) : V(_V), G(_G), C(_C), M(_M){
	match_prob[0] = 1.0f; match_prob[1] = M.ed/M.ee;

	unsigned max_nstates = 0;
	for (int h = 0 ; h < C.n_haplotypes ; h ++) if (C.indexes_pbwt_neighbour[h].size() > max_nstates) max_nstates = C.indexes_pbwt_neighbour[h].size();
	alpha = vector < aligned_vector32 < float > > (C.n_scaffold_variants, aligned_vector32 < float > (max_nstates, 0.0f));
	beta = aligned_vector32 < float > (max_nstates, 1.0f);
	Hvar.allocate(C.n_scaffold_variants, max_nstates);
	Hhap.allocate(max_nstates, C.n_scaffold_variants);
}

hmm_scaffold::~hmm_scaffold() {
	alpha.clear();
	beta.clear();
}

void hmm_scaffold::setup(unsigned int _hap) {
	hap = _hap;
	nstates = C.indexes_pbwt_neighbour[hap].size();

	Hvar.reallocateFast(C.n_scaffold_variants, nstates);
	Hhap.reallocateFast(nstates, C.n_scaffold_variants);

	Hhap.subset(C.Hhap, C.indexes_pbwt_neighbour[hap]);
	Hhap.transpose(Hvar);
}

double hmm_scaffold::forward() {
	float sum;
	double loglik = 0.0;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;
	const __m256i _vshift_count = _mm256_set_epi32(31,30,29,28,27,26,25,24);
	for (int vs = 0 ; vs < C.n_scaffold_variants ; vs ++) {
		const std::array<float,2> emit = {match_prob[C.Hhap.get(hap, vs)], match_prob[1-C.Hhap.get(hap, vs)]};
		const __m256 _emit0 = _mm256_set1_ps(emit[0]);
		const __m256 _emit1 = _mm256_set1_ps(emit[1]);

		if (!vs) {
			const float f0 = 1.0f / nstates;
			const __m256 _f0 = _mm256_set1_ps(f0);
			__m256 _sum = _mm256_set1_ps(0.0f);
			int offset = 0;
			for (int k = 0 ; k < nstatesMD8 ; k += 8) {
				const __m256i _mask = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )Hvar.getByte(vs, k)), _vshift_count);
				const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask));
				const __m256 _prob_curr = _mm256_mul_ps(_emiss, _f0);
				_sum = _mm256_add_ps(_sum, _prob_curr);
				_mm256_store_ps(&alpha[vs][k], _prob_curr);
				offset += 8;
			}
			sum = (offset > 0)?horizontal_add(_sum):0.0f;
			for (; offset < nstates ; offset ++) {
				alpha[vs][offset] = f0 * emit[Hvar.get(vs, offset)];
				sum += alpha[vs][offset];
			}
			loglik += log(sum);
		} else {
			const float f0 = M.t[vs-1] / nstates;
			const float f1 = M.nt[vs-1] / sum;
			const __m256 _f0 = _mm256_set1_ps(f0);
			const __m256 _f1 = _mm256_set1_ps(f1);
			__m256 _sum = _mm256_set1_ps(0.0f);
			int offset = 0;
			for (int k = 0 ; k < nstatesMD8 ; k += 8) {
				const __m256i _mask = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )Hvar.getByte(vs, k)), _vshift_count);
				const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask));
				const __m256 _prob_prev = _mm256_load_ps(&alpha[vs-1][k]);
				const __m256 _prob_temp = _mm256_fmadd_ps(_prob_prev, _f1, _f0);
				const __m256 _prob_curr = _mm256_mul_ps(_prob_temp, _emiss);
				_sum = _mm256_add_ps(_sum, _prob_curr);
				_mm256_store_ps(&alpha[vs][k], _prob_curr);
				offset += 8;
			}
			sum = (offset > 0)?horizontal_add(_sum):0.0f;
			for (; offset < nstates ; offset ++) {
				alpha[vs][offset] = (alpha[vs-1][offset]*f1+f0)*emit[Hvar.get(vs, offset)];
				sum += alpha[vs][offset];
			}
			loglik += log(sum);
		}
	}
	return loglik;
}

void hmm_scaffold::backward(vector < vector < unsigned int > > & cevents, vector < int > & vpath) {
	float sum = 0.0f, scale = 0.0f;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;
	const __m256i _vshift_count = _mm256_set_epi32(31,30,29,28,27,26,25,24);
	aligned_vector32 < float > alphaXbeta_curr = aligned_vector32 < float >(nstates, 0.0f);
	aligned_vector32 < float > alphaXbeta_prev = aligned_vector32 < float >(nstates, 0.0f);

	//vpath = vector < int > (C.n_scaffold_variants, -1);

	for (int vs = C.n_scaffold_variants - 1 ; vs >= 0 ; vs --) {

		//
		const std::array<float,2> emit = {match_prob[C.Hhap.get(hap, vs)], match_prob[1-C.Hhap.get(hap, vs)]};
		const __m256 _emit0 = _mm256_set1_ps(emit[0]);
		const __m256 _emit1 = _mm256_set1_ps(emit[1]);

		//Transitions
		if (vs == C.n_scaffold_variants - 1) fill (beta.begin(), beta.end(), 1.0f / nstates);
		else {
			const float f0 = M.t[vs] / nstates;
			const float f1 = M.nt[vs] / sum;
			const __m256 _f0 = _mm256_set1_ps(f0);
			const __m256 _f1 = _mm256_set1_ps(f1);
			int offset = 0;
			for (int k = 0 ; k < nstatesMD8 ; k += 8) {
				const __m256 _prob_prev = _mm256_load_ps(&beta[k]);
				const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _f1, _f0);
				_mm256_store_ps(&beta[k], _prob_curr);
				offset += 8;
			}
			for (; offset < nstates ; offset ++) beta[offset] = (beta[offset]*f1+f0);
		}

		//Products
		__m256 _scale = _mm256_set1_ps(0.0f);
		int offset = 0;
		for (int k = 0 ; k < nstatesMD8 ; k += 8) {
			const __m256 _prob_temp = _mm256_mul_ps(_mm256_load_ps(&alpha[vs][k]), _mm256_load_ps(&beta[k]));
			_mm256_store_ps(&alphaXbeta_curr[k], _prob_temp);
			_scale = _mm256_add_ps(_scale, _prob_temp);
			offset += 8;
		}
		scale = (offset > 0)?horizontal_add(_scale):0.0f;
		for (; offset < nstates ; offset ++) {
			alphaXbeta_curr[offset] = alpha[vs][offset] * beta[offset];
			scale += alphaXbeta_curr[offset];
		}
		scale = 1.0f / scale;
		_scale = _mm256_set1_ps(scale);
		offset = 0;
		for (int k = 0 ; k < nstatesMD8 ; k += 8) {
			const __m256 _prob_temp = _mm256_mul_ps(_mm256_load_ps(&alphaXbeta_curr[k]), _scale);
			_mm256_store_ps(&alphaXbeta_curr[k], _prob_temp);
			offset += 8;
		}
		for (; offset < nstates ; offset ++) alphaXbeta_curr[offset] *= scale;

		//Emission
		__m256 _sum = _mm256_set1_ps(0.0f);
		offset = 0;
		for (int k = 0 ; k < nstatesMD8 ; k += 8) {
			const __m256i _mask = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )Hvar.getByte(vs, k)), _vshift_count);
			const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask));
			const __m256 _prob_prev = _mm256_load_ps(&beta[k]);
			const __m256 _prob_curr = _mm256_mul_ps(_prob_prev, _emiss);
			_sum = _mm256_add_ps(_sum, _prob_curr);
			_mm256_store_ps(&beta[k], _prob_curr);
			offset += 8;
		}
		sum = (offset > 0)?horizontal_add(_sum):0.0f;
		for (; offset < nstates ; offset ++) {
			beta[offset] *= emit[Hvar.get(vs, offset)];
			sum += beta[offset];
		}

		//vpath[vs] = std::distance(alphaXbeta_curr.begin(),std::max_element(alphaXbeta_curr.begin(), alphaXbeta_curr.end()));

		//Storage
		if (cevents[vs+1].size()) {
			if (vs == C.n_scaffold_variants-1) copy(alphaXbeta_curr.begin(), alphaXbeta_curr.begin() + nstates, alphaXbeta_prev.begin());

			//Impute from full conditioning set
			for (int vr = 0 ; vr < cevents[vs+1].size() ; vr ++) {
				G.phaseLiAndStephens(cevents[vs+1][vr], hap, alphaXbeta_prev, alphaXbeta_curr, C.indexes_pbwt_neighbour[hap], 0.5001f);
			}
		}

		//Saving products
		copy(alphaXbeta_curr.begin(), alphaXbeta_curr.begin() + nstates, alphaXbeta_prev.begin());
	}

	if (cevents[0].size()) {

		//Impute from full conditioning set
		for (int vr = 0 ; vr < cevents[0].size() ; vr ++) {
			G.phaseLiAndStephens(cevents[0][vr], hap, alphaXbeta_curr, alphaXbeta_curr, C.indexes_pbwt_neighbour[hap], 0.5001f);
		}
	}
}

void hmm_scaffold::viterbi(vector < int > & path) {
	float sum, scale, maxv_prev, maxv_curr;
	int maxi_curr, maxi_prev;
	vector < vector < int > > _viterbi_paths = vector < vector < int > > (C.n_scaffold_variants, vector < int > (nstates, 0));
	vector < float > _viterbi_probs = vector < float > (nstates, 0.0f);

	//FORWARD PASS
	for (int vs = 0 ; vs < C.n_scaffold_variants ; vs ++) {
		const std::array < float, 2 > emit = { match_prob[C.Hhap.get(hap, vs)], match_prob[1-C.Hhap.get(hap, vs)] };

		if (!vs) {
			maxi_curr = -1;
			sum = maxv_curr = 0.0f;
			for (int k = 0; k < nstates ; k ++) {
				_viterbi_probs[k] = emit[Hvar.get(vs, k)];
				if (_viterbi_probs[k] > maxv_curr) {
					maxv_curr = _viterbi_probs[k];
					maxi_curr = k;
				}
				sum += _viterbi_probs[k];
			}
		} else {
			maxi_curr = -1;
			scale = 1.0f / sum;
			sum = maxv_curr = 0.0f;

			for (int k = 0 ; k < nstates ; k ++) {

				float prob_yrecomb = M.t[vs-1] * maxv_prev * scale;
				float prob_nrecomb = M.nt[vs-1] * _viterbi_probs[k] * scale;

				if (prob_yrecomb > prob_nrecomb) {		// I switch copying from the most likely state
					_viterbi_probs[k] = prob_yrecomb;
					_viterbi_paths[vs][k] = maxi_prev;
				} else {								// I stay copying from the same state
					_viterbi_probs[k] = prob_nrecomb;
					_viterbi_paths[vs][k] = k;
				}

				_viterbi_probs[k] *= emit[Hvar.get(vs, k)];

				if (_viterbi_probs[k] > maxv_curr) {
					maxv_curr = _viterbi_probs[k];
					maxi_curr = k;
				}

				sum += _viterbi_probs[k];
			}
		}

		maxi_prev = maxi_curr;
		maxv_prev = maxv_curr;
	}

	//BACKTRACKING PASS
	path = vector < int > (C.n_scaffold_variants, maxi_curr);
	for (int vs = C.n_scaffold_variants - 1 ; vs > 0; vs --)
		path[vs-1] = _viterbi_paths[vs][path[vs]];
}
