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

float rare_genotype::ee = 0.9999f;
float rare_genotype::ed = 0.0001f;

genotype_set::genotype_set() {
	clear();
}

genotype_set::~genotype_set() {
	clear();
}

void genotype_set::clear() {
	n_rare_variants = 0.0;
	n_samples = 0.0;
	names.clear();
	GRvar_genotypes.clear();
	GRind_genotypes.clear();
	nmiss_total = 0;
	nmiss_imputation = 0;
	nmiss_singleton = 0;
	nhets_total = 0;
	nhets_families = 0;
	nhets_imputation = 0;
	nhets_coalescent = 0;
}

void genotype_set::imputeMonomorphic() {
	for (int vr = 0 ; vr < n_rare_variants ; vr ++) {
		bool mono = true;
		for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {
			if (!GRvar_genotypes[vr][r].mis) mono = false;
		}
		if (mono) {
			nmiss_singleton += GRvar_genotypes[vr].size();
			GRvar_genotypes[vr].clear();
		}
	}
	vrb.bullet(stb.str(nmiss_singleton) + " missing genotypes imputed at monomorphic sites");
}

void genotype_set::allocate(variant_map & V, unsigned int _n_samples, unsigned int _n_scaffold_variants, unsigned int _n_rare_variants) {
	tac.clock();
	assert(_n_rare_variants > 0);

	n_scaffold_variants = _n_scaffold_variants;
	n_rare_variants = _n_rare_variants;
	n_samples = _n_samples;

	GRvar_genotypes = vector < vector < rare_genotype > > (n_rare_variants);
	GRind_genotypes = vector < vector < rare_genotype > > (n_samples);
	MAP_R2S = vector < unsigned int > (n_rare_variants);
	major_alleles = vector < bool > (n_rare_variants, false);
	for (int r = 0 ; r < V.sizeRare() ; r ++) major_alleles[r] = !V.vec_rare[r]->minor;

	//Mapping
	for (int vt = 0 ; vt < V.vec_scaffold[0]->idx_full ; vt ++) {
		if (V.vec_full[vt]->type == VARTYPE_RARE) MAP_R2S[V.vec_full[vt]->idx_rare] = 0;
	}
	for (int vs = 1 ; vs < V.sizeScaffold() ; vs ++) {
		for (int vt = V.vec_scaffold[vs-1]->idx_full+1 ; vt < V.vec_scaffold[vs]->idx_full ; vt ++) {
			if (V.vec_full[vt]->type == VARTYPE_RARE) MAP_R2S[V.vec_full[vt]->idx_rare] = vs;
		}
	}
	for (int vt = V.vec_scaffold.back()->idx_full + 1 ; vt < V.sizeFull() ; vt ++) {
		if (V.vec_full[vt]->type == VARTYPE_RARE) MAP_R2S[V.vec_full[vt]->idx_rare] = V.sizeScaffold();
	}

	vrb.bullet("Genotype set allocation [#rare=" + stb.str(n_rare_variants) + " / #samples=" + stb.str(n_samples) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_set::fillup_by_transpose_V2I() {
	tac.clock();

	//Clear-up
	for (int i = 0 ; i < n_samples ; i ++) GRind_genotypes[i].clear();

	//Fill-up
	for (int vr = 0 ; vr < n_rare_variants ; vr ++) {
		for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {
			unsigned int sample_idx = GRvar_genotypes[vr][r].idx;
			GRind_genotypes[sample_idx].push_back(GRvar_genotypes[vr][r]);
			GRind_genotypes[sample_idx].back().idx = vr;
		}
	}
	vrb.bullet("Genotype set transpose V2I (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_set::merge_by_transpose_I2V() {
	tac.clock();

	//Merge
	unsigned int ndone = 0;
	for (int i = 0 ; i < n_samples ; i ++) {
		for (int r = 0 ; r < GRind_genotypes[i].size() ; r ++) {
			unsigned int var_idx = GRind_genotypes[i][r].idx;
			bool found = false;
			for (int e = 0 ; e < GRvar_genotypes[var_idx].size() && !found; e ++ ) {
				if (GRvar_genotypes[var_idx][e].idx == i) {
					if (!GRvar_genotypes[var_idx][e].pha) {
						GRvar_genotypes[var_idx][e].al0 = GRind_genotypes[i][r].al0;
						GRvar_genotypes[var_idx][e].al1 = GRind_genotypes[i][r].al1;
						GRvar_genotypes[var_idx][e].prob = GRind_genotypes[i][r].prob;
						ndone++;
					}
					found = true;
				}
			}

			assert(found);
		}
	}
	vrb.bullet("Genotype set transpose merge I2V [n=" + stb.str(ndone) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_set::mapUnphasedOntoScaffold(int ind, vector < vector < unsigned int > > & map) {
	map.clear();
	map = vector < vector < unsigned int > > (n_scaffold_variants+1, vector < unsigned int > ());

	//Rare
	for (int vr = 0 ; vr < GRind_genotypes[ind].size() ; vr ++) {
		unsigned int idx = GRind_genotypes[ind][vr].idx;
		if (!GRind_genotypes[ind][vr].pha) map[MAP_R2S[idx]].push_back(idx);
	}
}


