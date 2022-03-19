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
	clear();
}

genotype_set::~genotype_set() {
	clear();
}

void genotype_set::clear() {
	n_rare_variants = 0.0;
	n_common_variants = 0.0;
	n_samples = 0.0;
	GCalleles.clear();
	GCmissing.clear();
	GRindexes.clear();
	GRhets.clear();
	GRmissing.clear();
}

void genotype_set::allocate(unsigned int _n_samples, unsigned int _n_rare_variants, unsigned int _n_common_variants) {
	tac.clock();

	n_rare_variants = _n_rare_variants;
	n_common_variants = _n_common_variants;
	n_samples = _n_samples;

	if (n_common_variants > 0) {
		GCalleles = vector < vector < bool > > (n_common_variants, vector < bool > (2 * n_samples, false));
		GCmissing = vector < vector < bool > > (n_common_variants, vector < bool > (n_samples, false));
	}

	if (n_rare_variants > 0) {
		GRindexes = vector < vector < unsigned int > > (n_rare_variants);
		GRhets = vector < vector < bool > > (n_rare_variants);
		GRmissing = vector < vector < bool > > (n_rare_variants);
		GRalleles = vector < vector < bool > > (n_rare_variants);
	}

	vrb.bullet("GEN allocation [#common=" + stb.str(n_common_variants) + " / #rare=" + stb.str(n_rare_variants) + " / #samples=" + stb.str(n_samples) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_set::getUnphasedIndexes(vector < unsigned int > & VC, vector < unsigned int > & VR, vector < unsigned int > & IDX) {
	IDX.clear();
	for (int vc = 0 ; vc < VC.size() ; vc ++) {
		for (int i = 0 ; i < n_samples ; i++) {
			if (GCmissing[VC[vc]][i] || GCalleles[VC[vc]][2*i+0] != GCalleles[VC[vc]][2*i+1]) IDX.push_back(i);
		}
	}
	for (int vr = 0 ; vr < VR.size() ; vr ++) {
		for (int i = 0 ; i < GRindexes[VR[vr]].size() ; i++) {
			if (GRhets[VR[vr]][i] || GRmissing[VR[vr]][i]) IDX.push_back(GRindexes[VR[vr]][i]);
		}
	}
	sort(IDX.begin(), IDX.end());
	IDX.erase(unique(IDX.begin(), IDX.end()), IDX.end());
}
