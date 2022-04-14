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
#include <models/pileup_caller/pileup_caller_header.h>

pileup_caller::pileup_caller(haplotype_set & _H, genotype_set & _G, variant_map & _V, int _min_mapQ, int _min_baseQ) : H(_H), G(_G), V(_V) {
	min_mapQ = _min_mapQ;
	min_baseQ = _min_baseQ;
	fai_fname = "";
	n_rhets_total = 0;
	n_rhets_pired = 0;
	n_bases_mismatch = 0;
	n_bases_total = 0;
	n_bases_lowqual = 0;
	n_bases_indel = 0;
	n_bases_match = 0;
	n_pirs_mismatch = 0;
	n_pirs_total = 0;
}

pileup_caller::~pileup_caller() {
	R.clear();
	if (fai) fai_destroy(fai);
}

void pileup_caller::loadFASTA(string ffasta) {
	setenv("REF_CACHE", "", 0);
	setenv("REF_PATH", "fake_value_so_no_download", 0);
	fai_fname = ffasta;
	fai = fai_load(fai_fname.c_str());
	if (fai == NULL) vrb.error("Error while loading fasta index file");
}

het pileup_caller::getScafHetsLeft(int ind, int vr) {
	//Backward ...
	het sg;
	for (int vt = V.vec_rare[G.GRind_genotypes[ind][vr].idx]->idx_full - 1; vt >= 0 ; vt --) {
		// ...until we hit a scaffold SNP
		if (V.vec_full[vt]->type == VARTYPE_SCAF && V.vec_full[vt]->isSNP()) {
			//Is it an het too?
			bool a0 = H.Hvar.get(V.vec_full[vt]->idx_scaffold, 2*ind+0);
			bool a1 = H.Hvar.get(V.vec_full[vt]->idx_scaffold, 2*ind+1);
			if (a0 != a1) {
				sg.idx = V.vec_full[vt]->idx_scaffold;
				sg.pos = V.vec_full[vt]->bp - 1;
				sg.a0 = a0;
				sg.a1 = a1;
				sg.ref = V.vec_full[vt]->ref[0];
				sg.alt = V.vec_full[vt]->alt[0];
				break;
			}
		}
	}
	return sg;
}

het pileup_caller::getScafHetsRight(int ind, int vr) {
	//Forward ...
	het sg;
	for (int vt = V.vec_rare[G.GRind_genotypes[ind][vr].idx]->idx_full + 1 ; vt < V.sizeFull() ; vt ++) {
		// ...until we hit a scaffold SNP
		if (V.vec_full[vt]->type == VARTYPE_SCAF && V.vec_full[vt]->isSNP()) {
			//Is it an het too?
			bool a0 = H.Hvar.get(V.vec_full[vt]->idx_scaffold, 2*ind+0);
			bool a1 = H.Hvar.get(V.vec_full[vt]->idx_scaffold, 2*ind+1);
			if (a0 != a1) {
				sg.idx = V.vec_full[vt]->idx_scaffold;
				sg.pos = V.vec_full[vt]->bp - 1;
				sg.a0 = a0;
				sg.a1 = a1;
				sg.ref = V.vec_full[vt]->ref[0];
				sg.alt = V.vec_full[vt]->alt[0];
				break;
			}
		}
	}
	return sg;
}

bool pileup_caller::phaseWithPIRs(int ind, int vr, het &lhet, het &thet, het &rhet) {
	if (R.size() == 0) return false;

	int support[2][2] = {};
	for (map < string, pir >::iterator itR = R.begin() ; itR != R.end() ; itR ++) {
		if (itR->second.t_obs && itR->second.l_obs && itR->second.r_obs) {
			if (itR->second.l_all == lhet.a0 && itR->second.r_all == rhet.a0) support[0][itR->second.t_all] ++;
			if (itR->second.l_all == lhet.a1 && itR->second.r_all == rhet.a1) support[1][itR->second.t_all] ++;
		} else if (itR->second.t_obs && itR->second.l_obs) {
			if (itR->second.l_all == lhet.a0) support[0][itR->second.t_all] ++;
			if (itR->second.l_all == lhet.a1) support[1][itR->second.t_all] ++;
		} else if (itR->second.t_obs && itR->second.r_obs) {
			if (itR->second.r_all == rhet.a0) support[0][itR->second.t_all] ++;
			if (itR->second.r_all == rhet.a1) support[1][itR->second.t_all] ++;
		}
	}

	float sum = support[0][0] + support[1][1] + support[0][1] + support[1][0];
	float prob_01 = (support[0][0] + support[1][1]) * 1.0f / sum;
	float prob_10 = (support[0][1] + support[1][0]) * 1.0f / sum;

	n_pirs_mismatch += min(support[0][0] + support[1][1], support[0][1] + support[1][0]);
	n_pirs_total += support[0][0] + support[1][1] + support[0][1] + support[1][0];

	if (prob_01 > 0.99) {
		thet.a0 = 0;
		thet.a1 = 1;
		return true;
	} else if (prob_10 > 0.99) {
		thet.a0 = 1;
		thet.a1 = 0;
		return true;
	}
	return false;
}
