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

#include <containers/genotype_set/genotype_set_header.h>

using namespace std;

void genotype_set::phaseLiAndStephens(unsigned int vr, unsigned int hap, aligned_vector32 < float > & alphaXbeta_prev, aligned_vector32 < float > & alphaXbeta_curr, vector < unsigned int > & H, float threshold) {
	float p[2] = { 0.0f };

	int tidx = -1;
	for (int k = 0, r = 0; (k<H.size()) || (r<GRvar_genotypes[vr].size()) ;) {
		int tind = (r<GRvar_genotypes[vr].size())?GRvar_genotypes[vr][r].idx:-1;
		int tmis = (r<GRvar_genotypes[vr].size())?GRvar_genotypes[vr][r].mis:-1;
		int cind = (k<H.size())?H[k]/2:-1;

		if (tind == hap/2) {
			tidx = r;
			r++;
		} else if (tind < 0) {
			p[major_alleles[vr]] += alphaXbeta_prev[k] * 0.5f + alphaXbeta_curr[k] * 0.5f;
			k++;
		} else if (cind < 0) {
			r++;
		} else if (tind < cind) {
			r++;
		} else if (tind==cind) {
			if (!tmis) p[!major_alleles[vr]] += alphaXbeta_prev[k] * 0.5f + alphaXbeta_curr[k] * 0.5f;
			k++;
		} else {
			p[major_alleles[vr]] += alphaXbeta_prev[k] * 0.5f + alphaXbeta_curr[k] * 0.5f;
			k++;
		}
	}

	assert(tidx>=0);
	assert((p[0]+p[1])>=0);

	if (!GRvar_genotypes[vr][tidx].pha) {
		if (hap%2 == 0) {
			assert(GRvar_genotypes[vr][tidx].prob < 0.0f);
			GRvar_genotypes[vr][tidx].prob = p[1] / (p[0] + p[1]);
		} else {
			assert(GRvar_genotypes[vr][tidx].prob >= 0.0f);
			if (haploids[GRvar_genotypes[vr][tidx].idx]) GRvar_genotypes[vr][tidx].impute(GRvar_genotypes[vr][tidx].prob, p[1] / (p[0] + p[1]));
			else GRvar_genotypes[vr][tidx].phase(GRvar_genotypes[vr][tidx].prob, p[1] / (p[0] + p[1]));
			assert(!isnan(GRvar_genotypes[vr][tidx].prob));
			assert(!isinf(GRvar_genotypes[vr][tidx].prob));
			if (GRvar_genotypes[vr][tidx].het && GRvar_genotypes[vr][tidx].prob < threshold) {
				GRvar_genotypes[vr][tidx].prob = -1.0f;
				GRvar_genotypes[vr][tidx].pha = 0;
			} else GRvar_genotypes[vr][tidx].pha = 1;
		}
	}
}

void genotype_set::phaseCoalescentViterbi(unsigned int ind, vector < int > & pathH0, vector < int > & pathH1, hmm_parameters & M) {
	 //
	 vector < int > starts0, ends0, starts1, ends1;
	 starts0.push_back(0);
	 starts1.push_back(0);
	 for (int l = 1 ; l < pathH0.size() ; l ++) {
		 if (pathH0[l-1] != pathH0[l]) {
			 ends0.push_back(l-1);
			 starts0.push_back(l);
		 }
		 if (pathH1[l-1] != pathH1[l]) {
			 ends1.push_back(l-1);
			 starts1.push_back(l);
		 }
	 }
	 ends0.push_back(pathH0.size() - 1);
	 ends1.push_back(pathH1.size() - 1);

	 //
	 vector < float > pathM0 = vector < float > (pathH0.size(), 0.0f);
	 vector < float > pathM1 = vector < float > (pathH1.size(), 0.0f);
	 for (int e = 0 ; e < starts0.size() ; e ++) {
		 float lengthCM = M.cm[ends0[e]] - M.cm[starts0[e]];
		 for (int l = starts0[e] ; l <= ends0[e] ; l++) pathM0[l] = lengthCM;
	 }
	 for (int e = 0 ; e < starts1.size() ; e ++) {
		 float lengthCM = M.cm[ends1[e]] - M.cm[starts1[e]];
		 for (int l = starts1[e] ; l <= ends1[e] ; l++) pathM1[l] = lengthCM;
	 }

	 //
	for (int vr = 0 ; vr < GRind_genotypes[ind].size() ; vr ++) {
		unsigned int idx = GRind_genotypes[ind][vr].idx;
		if (!GRind_genotypes[ind][vr].pha) {
			int index = MAP_R2S[idx];

			float w0, w1;
			if (index == 0) {
				w0 = pathM0[0];
				w1 = pathM1[0];
			} else if (index == pathH0.size()) {
				w0 = pathM0.back();
				w1 = pathM1.back();
			} else {
				w0 = max(pathM0[index-1], pathM0[index]);
				w1 = max(pathM1[index-1], pathM1[index]);
			}

			if (w0 > w1) {
				GRind_genotypes[ind][vr].al0 = 0;
				GRind_genotypes[ind][vr].al1 = 1;
			} else {
				GRind_genotypes[ind][vr].al0 = 1;
				GRind_genotypes[ind][vr].al1 = 0;
			}

			//GRind_genotypes[ind][vr].prob = max(w0, w1) / (w0+w1);
			GRind_genotypes[ind][vr].prob = 0.5f;
		}
	}
}
