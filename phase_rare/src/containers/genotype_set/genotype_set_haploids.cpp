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

void genotype_set::mapHaploidsAndResetHets(string fhap) {
	//DATA

	tac.clock();
	vrb.title("Mapping Haploid samples");

	//Reading file
	uint32_t nlines = 0;
	string buffer;
	map < string, int > samples_str2idx;
	map < string, int > :: iterator itS;
	for (int i = 0 ; i < n_samples ; i ++) samples_str2idx.insert(pair < string, int > (names[i], i));
	haploids = std::vector < bool > (n_samples, false);
	input_file fd_hap(fhap);
	if (fd_hap.fail()) vrb.error("Cannot open pedigree file");
	while (getline(fd_hap, buffer)) {
		itS = samples_str2idx.find(buffer);
		if (itS != samples_str2idx.end()) haploids[itS->second] = true;
		nlines ++;
	}
	fd_hap.close();
	vrb.bullet("PED file parsing [n=" + stb.str(nlines) + "]");

	//Counting haploids / diploids
	uint32_t n_diploids = 0, n_haploids = 0;
	for (int i = 0 ; i < n_samples ; i ++) {
		n_diploids += !haploids[i];
		n_haploids += haploids[i];
	}
	vrb.bullet("#haploids=" + stb.str(n_haploids) + " / #diploids=" + stb.str(n_diploids));

	//Zero out hets in Haploid samples
	uint64_t n_het = 0, n_hom = 0;
	for (int vr = 0 ; vr < n_rare_variants ; vr ++) {
		for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {
			unsigned int idx = GRvar_genotypes[vr][r].idx;
			if (haploids[idx]) {
				if (GRvar_genotypes[vr][r].het) {
					GRvar_genotypes[vr][r].het = 0;
					GRvar_genotypes[vr][r].mis = 1;
					n_het ++;
				} else if (!GRvar_genotypes[vr][r].mis) {
					n_hom ++;
				}
			}
		}
	}
	vrb.bullet("#reset_hets=" + stb.str(n_het) + " (" + stb.str(n_het * 100.0f / (n_het+n_hom), 3) + "% of minor allele genotypes)");
}

