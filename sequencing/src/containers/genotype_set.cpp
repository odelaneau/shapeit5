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
#include <containers/genotype_set.h>

genotype_set::genotype_set() {
	n_target_samples = 0;
	n_common_variants = 0;
}

genotype_set::~genotype_set() {
	for (int i = 0 ; i< vecG.size() ; i ++) delete vecG[i];
	vecG.clear();
	n_target_samples = 0;
	n_common_variants = 0;
}

void genotype_set::allocate(unsigned int _n_target_samples, unsigned int _n_common_variants) {
	n_target_samples = _n_target_samples;
	n_common_variants = _n_common_variants;
	vecG = vector < genotype * > (n_target_samples);
	for (unsigned int i = 0 ; i < n_target_samples ; i ++) vecG[i] = new genotype (i, n_common_variants);
}


void genotype_set::imputeMonomorphic(variant_map & V) {
	for (unsigned int vt = 0, vc = 0 ; vt < V.sizeFull() ; vt ++) {
		if (!V.vec_full[vt]->rare) {
			if (V.vec_full[vt]->isMonomorphic()) {
				bool uallele = (V.vec_full[vt]->cref)?false:true;
				for (unsigned int i = 0 ; i < vecG.size() ; i ++) {
					VAR_SET_HOM(MOD2(vc), vecG[i]->Variants[DIV2(vc)]);
					uallele?VAR_SET_HAP0(MOD2(vc), vecG[i]->Variants[DIV2(vc)]):VAR_CLR_HAP0(MOD2(vc), vecG[i]->Variants[DIV2(vc)]);
					uallele?VAR_SET_HAP1(MOD2(vc), vecG[i]->Variants[DIV2(vc)]):VAR_CLR_HAP1(MOD2(vc), vecG[i]->Variants[DIV2(vc)]);
				}
				if (uallele) V.vec_full[vt]->cref = 0;
				else V.vec_full[vt]->calt = 0;
				V.vec_full[vt]->cmis = 0;
			}
			vc++;
		} else {
			if (V.vec_full[vt]->isMonomorphic() && V.vec_full[vt]->cmis > 0) {
				bool uallele = (V.vec_full[vt]->cref)?false:true;
				for (unsigned int i = 0 ; i < vecG.size() ; i ++) {
					int idx_rare = -1;
					for (int r = 0 ; r < vecG[i]->RareIndexes.size() && idx_rare < 0 ; r ++) if (vecG[i]->RareIndexes[r] == vt) idx_rare = r;
					if (idx_rare >= 0) {
						VAR_SET_HOM(MOD2(idx_rare), vecG[i]->RareVariants[DIV2(idx_rare)]);
						uallele?VAR_SET_HAP0(MOD2(idx_rare), vecG[i]->RareVariants[DIV2(idx_rare)]):VAR_CLR_HAP0(MOD2(idx_rare), vecG[i]->RareVariants[DIV2(idx_rare)]);
						uallele?VAR_SET_HAP1(MOD2(idx_rare), vecG[i]->RareVariants[DIV2(idx_rare)]):VAR_CLR_HAP1(MOD2(idx_rare), vecG[i]->RareVariants[DIV2(idx_rare)]);
					}
				}
				if (uallele) V.vec_full[vt]->cref = 0;
				else V.vec_full[vt]->calt = 0;
				V.vec_full[vt]->cmis = 0;
			}
		}
	}
}

unsigned int genotype_set::largestNumberOfTransitions() {
	unsigned int maxT = 0;
	for (int i = 0 ; i < n_target_samples ; i ++) {
		unsigned int nTrans = vecG[i]->n_transitions;
		if (nTrans > maxT) maxT = nTrans;
	}
	return maxT;
}

unsigned int genotype_set::largestNumberOfMissings() {
	unsigned int maxM = 0;
	for (int i = 0 ; i < n_target_samples ; i ++) {
		unsigned int nMis = vecG[i]->n_missing;
		if (nMis> maxM) maxM = nMis;
	}
	return maxM;
}

unsigned int genotype_set::largestNumberOfRares() {
	unsigned int maxR = 0;
	for (int i = 0 ; i < n_target_samples ; i ++) {
		unsigned int nRare = vecG[i]->RareIndexes.size();
		if (nRare> maxR) maxR = nRare;
	}
	return maxR;
}

unsigned long genotype_set::numberOfSegments() {
	unsigned long size = 0;
	for (int i = 0 ; i < n_target_samples ; i ++) size += vecG[i]->n_segments;
	return size;
}

void genotype_set::masking() {
	for (int i = 0 ; i < n_target_samples ; i ++) vecG[i]->mask();
}

void genotype_set::solve() {
	tac.clock();
	for (int i = 0 ; i < vecG.size() ; i ++) vecG[i]->solve();
	vrb.bullet("HAP solving (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
