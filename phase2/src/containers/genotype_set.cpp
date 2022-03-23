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
	names.clear();
	GRvar_genotypes.clear();
	GRind_genotypes.clear();
}

//Mapping on scaffold
vector < int > MAPC_vs_left;
vector < int > MAPC_vs_right;
vector < int > MAPR_vs_left;
vector < int > MAPR_vs_right;


void genotype_set::allocate(variant_map & V, unsigned int _n_samples, unsigned int _n_scaffold_variants, unsigned int _n_rare_variants, unsigned int _n_common_variants) {
	tac.clock();

	n_scaffold_variants = _n_scaffold_variants;
	n_rare_variants = _n_rare_variants;
	n_common_variants = _n_common_variants;
	n_samples = _n_samples;

	if (n_common_variants > 0) {
		GCvar_alleles.allocate(n_common_variants, 2 * n_samples);
		GCind_alleles.allocate(2 * n_samples, n_common_variants);
		GCvar_missing.allocate(n_common_variants, n_samples);
		GCind_missing.allocate(n_samples, n_common_variants);
		MAPC_vs_left = vector < unsigned int > (n_common_variants, -1);
		MAPC_vs_right = vector < unsigned int > (n_common_variants, -1);

		GCvar_truth.allocate(n_common_variants, 2 * n_samples);
	}

	if (n_rare_variants > 0) {
		GRvar_genotypes = vector < vector < rare_genotype > > (n_rare_variants);
		GRind_genotypes = vector < vector < rare_genotype > > (n_samples);
		MAPR_vs_left = vector < unsigned int > (n_rare_variants, -1);
		MAPR_vs_right = vector < unsigned int > (n_rare_variants, -1);
		major_alleles = vector < bool > (n_rare_variants, false);
		for (int r = 0 ; r < V.sizeRare() ; r ++) major_alleles[r] = !V.vec_rare[r]->minor;

		GRvar_truth = vector < vector < bool > > (n_rare_variants);

	}

	//Mapping
	for (int vt = 0 ; vt < V.vec_scaffold[0]->idx_full ; vt ++) {
		if (V.vec_full[vt]->type == VARTYPE_RARE) { MAPR_vs_left[V.vec_full[vt]->idx_rare] = 0; MAPR_vs_right[V.vec_full[vt]->idx_rare] = 1; }
		if (V.vec_full[vt]->type == VARTYPE_COMM) { MAPC_vs_left[V.vec_full[vt]->idx_common] = 0; MAPC_vs_right[V.vec_full[vt]->idx_common] = 1; }
	}
	for (int vs = 1 ; vs < V.sizeScaffold() ; vs ++) {
		for (int vt = V.vec_scaffold[vs-1]->idx_full+1 ; vt < V.vec_scaffold[vs]->idx_full ; vt ++) {
			if (V.vec_full[vt]->type == VARTYPE_RARE) { MAPR_vs_left[V.vec_full[vt]->idx_rare] = vs; MAPR_vs_right[V.vec_full[vt]->idx_rare] = vs+1; }
			if (V.vec_full[vt]->type == VARTYPE_COMM) { MAPC_vs_left[V.vec_full[vt]->idx_common] = vs; MAPC_vs_right[V.vec_full[vt]->idx_common] = vs+1; }
		}
	}
	for (int vt = V.vec_scaffold.back()->idx_full ; vt < V.sizeFull() ; vt ++) {
		if (V.vec_full[vt]->type == VARTYPE_RARE) { MAPR_vs_left[V.vec_full[vt]->idx_rare] = V.sizeScaffold()-1; MAPR_vs_right[V.vec_full[vt]->idx_rare] = V.sizeScaffold(); }
		if (V.vec_full[vt]->type == VARTYPE_COMM) { MAPC_vs_left[V.vec_full[vt]->idx_common] = V.sizeScaffold()-1; MAPC_vs_right[V.vec_full[vt]->idx_common] = V.sizeScaffold(); }
	}

	/*
	for (int vr = 0 ; vr < n_rare_variants ; vr ++)
		cout << "R " << vr << " " << MAPR_vs_left[vr] << " " << MAPR_vs_right[vr] << endl;
	for (int vc = 0 ; vc < n_common_variants ; vc ++)
		cout << "C " << vc << " " << MAPC_vs_left[vc] << " " << MAPC_vs_right[vc] << endl;
	 */

	vrb.bullet("Genotype set allocation [#common=" + stb.str(n_common_variants) + " / #rare=" + stb.str(n_rare_variants) + " / #samples=" + stb.str(n_samples) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_set::transpose() {
	tac.clock();
	if (n_common_variants > 0) {
		GCvar_alleles.transpose(GCind_alleles, n_common_variants, 2 * n_samples);
		GCvar_missing.transpose(GCind_missing, n_common_variants, n_samples);
	}

	if (n_rare_variants > 0) {
		for (int vr = 0 ; vr < n_rare_variants ; vr ++) {
			for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {
				unsigned int sample_idx = GRvar_genotypes[vr][r].idx;
				GRind_genotypes[sample_idx].push_back(GRvar_genotypes[vr][r]);
				GRind_genotypes[sample_idx].back().idx = vr;
			}
		}
	}
	vrb.bullet("Genotype set transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_set::mapUnphasedOntoScaffold(int ind, vector < bool > & map) {
	map.clear();
	map = vector < bool > (n_scaffold_variants+1, false);

	//Common
	for (int vc = 0 ; vc < n_common_variants ; vc ++) {
		if (GCind_missing.get(ind, vc) || GCind_alleles.get(2*ind+0, vc) != GCind_alleles.get(2*ind+1, vc)) {
			map[MAPC_vs_left[vc]] = true;
		}
	}
	//Rare
	for (int vr = 0 ; vr < GRind_genotypes[ind].size() ; vr ++) {
		if (GRind_genotypes[ind][vr].mis || GRind_genotypes[ind][vr].het) {
			map[MAPR_vs_left[GRind_genotypes[ind][vr].idx]] = true;
		}
	}
}
