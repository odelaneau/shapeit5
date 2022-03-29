/*******************************************************************************
 * Copyright (C) 2018-2022 Olivier Delaneau, University of Lausanne
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
#ifndef _HAPLOTYPE_SEGMENT_SINGLE_H
#define _HAPLOTYPE_SEGMENT_SINGLE_H

#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>

#ifdef __AVX2__

#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

#endif

class haplotype_segment {
private:
	//EXTERNAL DATA
	hmm_parameters & M;
	genotype * G;
	bitmatrix Hhap, Hvar;

	//COORDINATES & CONSTANTS
	int segment_first;
	int segment_last;
	int locus_first;
	int locus_last;
	int ambiguous_first;
	int ambiguous_last;
	int missing_first;
	int missing_last;
	int transition_first;
	int transition_last;
	unsigned int n_cond_haps;
	unsigned int n_missing;

	//CURSORS
	int curr_segment_index;
	int curr_segment_locus;
	int curr_abs_locus;
	int prev_abs_locus;
	int curr_rel_locus;
	int curr_rel_locus_offset;
	int curr_abs_ambiguous;
	int curr_abs_transition;
	int curr_abs_missing;
	int curr_rel_missing;


	//DYNAMIC ARRAYS
	float probSumT;
#ifdef __AVX2__
	aligned_vector32 < float > prob;
	aligned_vector32 < float > probSumK;
	aligned_vector32 < float > probSumH;
	vector < aligned_vector32 < float > > Alpha;
	vector < aligned_vector32 < float > > AlphaSum;
	vector < int > AlphaLocus;
	aligned_vector32 < float > AlphaSumSum;
	vector < aligned_vector32 < float > > AlphaMissing;
	vector < aligned_vector32 < float > > AlphaSumMissing;
	float HProbs [HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));
	double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));
#else
	vector < float > prob;
	vector < float > probSumK;
	vector < float > probSumH;
	vector < vector < float > > Alpha;
	vector < vector < float > > AlphaSum;
	vector < int > AlphaLocus;
	vector < float > AlphaSumSum;
	vector < vector < float > > AlphaMissing;
	vector < vector < float > > AlphaSumMissing;
	float HProbs [HAP_NUMBER * HAP_NUMBER];
	double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER];
#endif

	//STATIC ARRAYS
	float sumHProbs;
	double sumDProbs;
	float g0[HAP_NUMBER], g1[HAP_NUMBER];
	float nt, yt;

	//INLINED AND UNROLLED ROUTINES
	void INIT_HOM();
	void INIT_AMB();
	void INIT_MIS();
	bool RUN_HOM(char);
	void RUN_AMB();
	void RUN_MIS();
	void COLLAPSE_HOM();
	void COLLAPSE_AMB();
	void COLLAPSE_MIS();
	void SUMK();
	void IMPUTE(vector < float > & );
	bool TRANS_HAP();
	bool TRANS_DIP_MULT();
	bool TRANS_DIP_ADD();
	void SET_FIRST_TRANS(vector < double > & );
	int SET_OTHER_TRANS(vector < double > & );

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_segment(genotype *, bitmatrix &, vector < unsigned int > &, window &, hmm_parameters &);
	~haplotype_segment();

	//void fetch();
	void forward();
	int backward(vector < double > &, vector < float > &);
};

/*******************************************************************************/
/*****************			HOMOZYGOUS GENOTYPE			************************/
/*******************************************************************************/

#ifdef __AVX2__
inline
void haplotype_segment::INIT_HOM() {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	__m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _mm256_set1_ps((ag==ah)?1.0f:M.ed/M.ee);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#else
inline
void haplotype_segment::INIT_HOM() {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (ag==ah)?1.0f:M.ed/M.ee);
		for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#endif

#ifdef __AVX2__
inline
bool haplotype_segment::RUN_HOM(char rare_allele) {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	if (rare_allele < 0 || ag == rare_allele) {
		__m256 _sum = _mm256_set1_ps(0.0f);
		//__m256 _factor = _mm256_set1_ps(M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
		__m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
		__m256 _tFreq = _mm256_load_ps(&probSumH[0]);
		_tFreq = _mm256_mul_ps(_tFreq, _factor);
		//__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT);
		__m256 _nt = _mm256_set1_ps(nt / probSumT);
		__m256 _mismatch = _mm256_set1_ps(M.ed/M.ee);
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
			__m256 _prob = _mm256_load_ps(&prob[i]);
			_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
			if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
			_sum = _mm256_add_ps(_sum, _prob);
			_mm256_store_ps(&prob[i], _prob);
		}
		_mm256_store_ps(&probSumH[0], _sum);
		probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
		return true;
	}
	return false;
}
#else
inline
bool haplotype_segment::RUN_HOM(char rare_allele) {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	if (rare_allele < 0 || ag == rare_allele) {
		//vector < float > _tFreq = vector < float >(HAP_NUMBER, M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
		vector < float > _tFreq = vector < float >(HAP_NUMBER, yt / (n_cond_haps * probSumT));
		for (int h = 0 ; h < HAP_NUMBER ; h++) _tFreq[h] *= probSumH[h];
		//float _nt = M.nt[curr_abs_locus-forward] / probSumT;
		float _nt = nt / probSumT;
		float _mismatch = M.ed/M.ee;
		fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
			for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] = (prob[i+h]*_nt)+_tFreq[h];
			if (ag!=ah) for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] *= _mismatch;
			for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
		}
		probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
		return true;
	}
	return false;
}
#endif

#ifdef __AVX2__
inline
void haplotype_segment::COLLAPSE_HOM() {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	__m256 _sum = _mm256_set1_ps(0.0f);
	//__m256 _tFreq = _mm256_set1_ps(M.t[curr_abs_locus-forward] / n_cond_haps);
	__m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);					////Check divide by probSumT here!
	//__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	__m256 _mismatch = _mm256_set1_ps(M.ed/M.ee);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _mm256_set1_ps(probSumK[k]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#else
inline
void haplotype_segment::COLLAPSE_HOM() {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	//float _tFreq = M.t[curr_abs_locus-forward] / n_cond_haps;
	float _tFreq = yt / n_cond_haps;					////Check divide by probSumT here!
	//float _nt = M.nt[curr_abs_locus-forward] / probSumT;
	float _nt = nt / probSumT;
	float _mismatch = M.ed/M.ee;
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (ag==ah)?((probSumK[k]*_nt)+_tFreq):(((probSumK[k]*_nt)+_tFreq)*_mismatch));
		for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#endif

/*******************************************************************************/
/*****************			HETEROZYGOUS GENOTYPE			********************/
/*******************************************************************************/

#ifdef __AVX2__
inline
void haplotype_segment::INIT_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	__m256 _sum = _mm256_set1_ps(0.0f);
	__m256 _emit0 = _mm256_loadu_ps(&g0[0]);
	__m256 _emit1 = _mm256_loadu_ps(&g1[0]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = ah?_emit1:_emit0;
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#else
inline
void haplotype_segment::INIT_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		if (ah) memcpy(&prob[i], &g1[0], HAP_NUMBER*sizeof(float));
		else memcpy(&prob[i], &g0[0], HAP_NUMBER*sizeof(float));
		for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#endif

#ifdef __AVX2__
inline
void haplotype_segment::RUN_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	__m256 _sum = _mm256_set1_ps(0.0f);
	//__m256 _factor = _mm256_set1_ps(M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
	__m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
	__m256 _tFreq = _mm256_load_ps(&probSumH[0]);
	_tFreq = _mm256_mul_ps(_tFreq, _factor);
	//__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	__m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _mm256_load_ps(&prob[i]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_prob = _mm256_mul_ps(_prob, _emit[ah]);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#else
inline
void haplotype_segment::RUN_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	//vector < float > _tFreq = vector < float >(HAP_NUMBER, M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
	vector < float > _tFreq = vector < float >(HAP_NUMBER, yt / (n_cond_haps * probSumT));
	for (int h = 0 ; h < HAP_NUMBER ; h++) _tFreq[h] *= probSumH[h];
	//float _nt = M.nt[curr_abs_locus-forward] / probSumT;
	float _nt = nt / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] = (prob[i+h]*_nt)+_tFreq[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] *= ah?g1[h]:g0[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#endif

#ifdef __AVX2__
inline
void haplotype_segment::COLLAPSE_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	__m256 _sum = _mm256_set1_ps(0.0f);
	//__m256 _tFreq = _mm256_set1_ps(M.t[curr_abs_locus-forward] / n_cond_haps);
	__m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
	//__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	__m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _mm256_set1_ps(probSumK[k]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_prob = _mm256_mul_ps(_prob, _emit[ah]);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#else
inline
void haplotype_segment::COLLAPSE_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	//float _tFreq = M.t[curr_abs_locus-forward] / n_cond_haps;
	float _tFreq = yt / n_cond_haps;
	//float _nt = M.nt[curr_abs_locus-forward] / probSumT;
	float _nt = nt / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (probSumK[k]*_nt)+_tFreq);
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] *= ah?g1[h]:g0[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#endif

/*******************************************************************************/
/*****************			MISSING GENOTYPE			************************/
/*******************************************************************************/

inline
void haplotype_segment::INIT_MIS() {
	fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * n_cond_haps));
	fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
	probSumT = 1.0f;
}

#ifdef __AVX2__
inline
void haplotype_segment::RUN_MIS() {
	__m256 _sum = _mm256_set1_ps(0.0f);
	//__m256 _factor = _mm256_set1_ps(M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
	__m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
	__m256 _tFreq = _mm256_load_ps(&probSumH[0]);
	_tFreq = _mm256_mul_ps(_tFreq, _factor);
	//__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		__m256 _prob = _mm256_load_ps(&prob[i]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#else
inline
void haplotype_segment::RUN_MIS() {
	//vector < float > _tFreq = vector < float >(HAP_NUMBER, M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
	vector < float > _tFreq = vector < float >(HAP_NUMBER, yt / (n_cond_haps * probSumT));
	for (int h = 0 ; h < HAP_NUMBER ; h++) _tFreq[h] *= probSumH[h];
	//float _nt = M.nt[curr_abs_locus-forward] / probSumT;
	float _nt = nt / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] = (prob[i+h]*_nt)+_tFreq[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#endif

#ifdef __AVX2__
inline
void haplotype_segment::COLLAPSE_MIS() {
	__m256 _sum = _mm256_set1_ps(0.0f);
	//__m256 _tFreq = _mm256_set1_ps(M.t[curr_abs_locus-forward] / n_cond_haps);
	__m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
	//__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		__m256 _prob = _mm256_set1_ps(probSumK[k]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#else
inline
void haplotype_segment::COLLAPSE_MIS() {
	//float _tFreq = M.t[curr_abs_locus-forward] / n_cond_haps;
	float _tFreq = yt / n_cond_haps;
	//float _nt = M.nt[curr_abs_locus-forward] / probSumT;
	float _nt = nt / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (probSumK[k]*_nt)+_tFreq);
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
#endif

/*******************************************************************************/
/*****************					SUM Ks				************************/
/*******************************************************************************/

inline
void haplotype_segment::SUMK() {
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		probSumK[k] = prob[i+0] + prob[i+1] + prob[i+2] + prob[i+3] + prob[i+4] + prob[i+5] + prob[i+6] + prob[i+7];
	}
}

/*******************************************************************************/
/*****************		TRANSITION COMPUTATIONS			************************/
/*******************************************************************************/

#ifdef __AVX2__
inline
bool haplotype_segment::TRANS_HAP() {
	sumHProbs = 0.0f;
	unsigned int  curr_rel_segment_index = curr_segment_index-segment_first;
	yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1], prev_abs_locus);
	nt = 1.0f - yt;
	//float fact1 = M.nt[curr_abs_locus-1] / AlphaSumSum[curr_rel_segment_index - 1];
	float fact1 = nt / AlphaSumSum[curr_rel_segment_index - 1];
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		__m256 _sum = _mm256_set1_ps(0.0f);
		//float fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * M.t[curr_abs_locus - 1] / n_cond_haps;
		float fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * yt / n_cond_haps;
		for (int k = 0 ; k < n_cond_haps ; k ++) {
			__m256 _alpha = _mm256_set1_ps(Alpha[curr_rel_segment_index-1][k*HAP_NUMBER + h1] * fact1 + fact2);
			__m256 _beta = _mm256_load_ps(&prob[k*HAP_NUMBER]);
			_sum = _mm256_add_ps(_sum, _mm256_mul_ps(_alpha, _beta));
		}
		_mm256_store_ps(&HProbs[h1*HAP_NUMBER], _sum);
		sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
	}
	return (isnan(sumHProbs) || isinf(sumHProbs) || sumHProbs < numeric_limits<float>::min());
}
#else
inline
bool haplotype_segment::TRANS_HAP() {
	sumHProbs = 0.0f;
	unsigned int  curr_rel_segment_index = curr_segment_index-segment_first;
	yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1], curr_abs_locus);
	nt = 1.0f - yt;
	//float fact1 = M.nt[curr_abs_locus-1] / AlphaSumSum[curr_rel_segment_index - 1];
	float fact1 = nt / AlphaSumSum[curr_rel_segment_index - 1];
	fill_n(HProbs, HAP_NUMBER*HAP_NUMBER, 0.0f);
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		//float fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * M.t[curr_abs_locus - 1] / n_cond_haps;
		float fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * yt / n_cond_haps;
		for (int k = 0 ; k < n_cond_haps ; k ++) {
			for (int h2 = 0 ; h2 < HAP_NUMBER ; h2++) HProbs[h1*HAP_NUMBER+h2]+=((Alpha[curr_rel_segment_index-1][k*HAP_NUMBER + h1]*fact1 + fact2)*prob[k*HAP_NUMBER+h2]);
		}
		sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
	}
	return (isnan(sumHProbs) || isinf(sumHProbs) || sumHProbs < numeric_limits<float>::min());
}
#endif

inline
bool haplotype_segment::TRANS_DIP_MULT() {
	sumDProbs= 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) * ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (isnan(sumDProbs) || isinf(sumDProbs) || sumDProbs < numeric_limits<double>::min());
}

inline
bool haplotype_segment::TRANS_DIP_ADD() {
	sumDProbs = 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) + ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (isnan(sumDProbs) || isinf(sumDProbs) || sumDProbs < numeric_limits<double>::min());
}

#ifdef __AVX2__
inline
void haplotype_segment::IMPUTE(vector < float > & missing_probabilities) {
	__m256 _sum = _mm256_set1_ps(0.0f);
	__m256 _sumA [2]; _sumA[0] = _mm256_set1_ps(0.0f); _sumA[1] = _mm256_set1_ps(0.0f);
	__m256 _alphaSum = _mm256_load_ps(&AlphaSumMissing[curr_rel_missing][0]);
	__m256 _ones = _mm256_set1_ps(1.0f);
	_alphaSum = _mm256_div_ps(_ones, _alphaSum);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _mm256_load_ps(&prob[i]);
		__m256 _alpha = _mm256_load_ps(&AlphaMissing[curr_rel_missing][i]);
		_sum = _mm256_mul_ps(_mm256_mul_ps(_alpha, _alphaSum), _prob);
		_sumA[ah] = _mm256_add_ps(_sumA[ah], _sum);
	}
	float * prob0 = (float*)&_sumA[0];
	float * prob1 = (float*)&_sumA[1];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = prob1[h] / (prob0[h]+prob1[h]);
	}
}
#else
inline
void haplotype_segment::IMPUTE(vector < float > & missing_probabilities) {
	vector < vector < float > > _sumA = vector < vector < float > > (2, vector < float > (HAP_NUMBER, 0.0f));
	vector < float > _scale = AlphaSumMissing[curr_rel_missing];
	for (int h = 0 ; h < HAP_NUMBER ; h++) _scale[h] = 1.0f / _scale[h];
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		for (int h = 0 ; h < HAP_NUMBER ; h++) _sumA[ah][h] += (AlphaMissing[curr_rel_missing][i+h]*_scale[h])*prob[i+h];
	}
	for (int h = 0 ; h < HAP_NUMBER ; h ++) missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = _sumA[1][h] / (_sumA[0][h]+_sumA[1][h]);
}
#endif

#endif
