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
#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	clear();
}

haplotype_set::~haplotype_set() {
	clear();
}

void haplotype_set::clear() {
	n_site = 0;
	n_hap = 0;
	n_ind = 0;
}

void haplotype_set::allocate(unsigned long n_main_samples, unsigned long n_ref_samples, unsigned long n_variants) {
	n_ind = n_main_samples;
	n_hap = 2 * (n_main_samples + n_ref_samples);
	n_site = n_variants;
	H_opt_var.allocate(n_site, n_hap);
	H_opt_hap.allocate(n_hap, n_site);
}

void haplotype_set::updateHaplotypes(genotype_set & G, bool first_time) {
	tac.clock();
	for (unsigned int i = 0 ; i < G.n_ind ; i ++) {
		for (unsigned int v = 0 ; v < n_site ; v ++) {
			if (first_time || (VAR_GET_HET(MOD2(v), G.vecG[i]->Variants[DIV2(v)])) || (VAR_GET_MIS(MOD2(v), G.vecG[i]->Variants[DIV2(v)]))) {
				bool a0 = VAR_GET_HAP0(MOD2(v), G.vecG[i]->Variants[DIV2(v)]);
				bool a1 = VAR_GET_HAP1(MOD2(v), G.vecG[i]->Variants[DIV2(v)]);
				H_opt_hap.set(2*i+0, v, a0);
				H_opt_hap.set(2*i+1, v, a1);
			}
		}
	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_H2V(bool full, bool verbose) {
	if (verbose) tac.clock();
	if (!full) H_opt_hap.transpose(H_opt_var, 2*n_ind, n_site);
	else H_opt_hap.transpose(H_opt_var, n_hap, n_site);
	if (verbose) vrb.bullet("H2V transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_V2H(bool full, bool verbose) {
	if (verbose) tac.clock();
	if (!full) H_opt_var.transpose(H_opt_hap, n_site, 2*n_ind);
	else H_opt_var.transpose(H_opt_hap, n_site, n_hap);
	if (verbose) vrb.bullet("V2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

