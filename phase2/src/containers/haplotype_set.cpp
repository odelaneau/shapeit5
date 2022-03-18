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
#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	clear();
}

haplotype_set::~haplotype_set() {
	clear();
}

void haplotype_set::clear() {
	n_scaffold_variants = 0.0;
	n_rare_variants = 0.0;
	n_common_variants = 0.0;
	n_total_variants = 0.0;
	n_haplotypes = 0.0;
	n_samples = 0.0;
}

void haplotype_set::allocate(unsigned int _n_samples, unsigned int _n_scaffold_variants, unsigned int _n_rare_variants, unsigned int _n_common_variants, variant_map &) {
	tac.clock();

	n_scaffold_variants = _n_scaffold_variants;
	n_rare_variants = _n_rare_variants;
	n_common_variants = _n_common_variants;
	n_total_variants = n_scaffold_variants + n_rare_variants + n_common_variants;
	n_haplotypes = 2 * _n_samples;
	n_samples = _n_samples;

	var_type = vector < char > (n_total_variants, -1);
	for (int v = 0 ; v < n_total_variants ; v ++) var_type[v] = V.vec_pos[v]->type;

	HSvar.allocate(n_scaffold_variants, n_haplotypes);
	HShap.allocate(n_haplotypes, n_scaffold_variants);

	if (n_common_variants > 0) {
		HCvar.allocate(n_common_variants, n_haplotypes);
		HChap.allocate(n_haplotypes, n_common_variants);
	}

	if (n_rare_variants > 0) {
		HRvar = vector < vector < unsigned int > > (n_rare_variants, vector < unsigned int > ());
		HRhap = vector < vector < unsigned int > > (n_haplotypes, vector < unsigned int > ());
	}

	vrb.bullet("HAP allocation [#scaffold=" + stb.str(n_scaffold_variants) + " / #common=" + stb.str(n_common_variants) + " / #rare=" + stb.str(n_rare_variants) + " / #samples=" + stb.str(n_samples) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_H2V() {
	tac.clock();

	HShap.transpose(HSvar, n_haplotypes, n_scaffold_variants);

	if (n_common_variants > 0) HChap.transpose(HCvar, n_haplotypes, n_common_variants);

	if (n_rare_variants > 0) {
		for (unsigned int vr = 0 ; vr < n_rare_variants ; vr ++) HRvar[vr].clear();
		for (unsigned int h = 0 ; h < n_haplotypes ; h ++) for (unsigned int r = 0 ; r < HRhap[h].size() ; r ++) HRvar[HRhap[h][r]].push_back(h);
	}

	vrb.bullet("H2V transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_V2H(bool full) {
	tac.clock();

	HSvar.transpose(HShap, n_scaffold_variants, n_haplotypes);

	if (n_common_variants > 0) HCvar.transpose(HChap, n_common_variants, n_haplotypes);

	if (n_rare_variants > 0) {
		for (unsigned int h = 0 ; h < n_haplotypes ; h ++) HRhap[h].clear();
		for (unsigned int vr = 0 ; vr < n_rare_variants ; vr ++) for (unsigned int r = 0 ; r < HRvar[vr].size() ; r ++) HRhap[HRvar[vr][r]].push_back(vr);
	}

	vrb.bullet("V2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

