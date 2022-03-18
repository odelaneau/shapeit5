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

#include <containers/conditioning_set/conditioning_set_header.h>

void conditioning_set::select() {
	tac.clock();
	vrb.progress("  * PBWT selection", 0.0f);

	vector < int > A = vector < int > (n_haplotypes, 0);
	vector < int > B = vector < int > (n_haplotypes, 0);
	vector < int > C = vector < int > (n_haplotypes, 0);
	vector < int > D = vector < int > (n_haplotypes, 0);
	vector < int > M = vector < int > (n_haplotypes, -1);
	iota(A.begin(), A.end(), 0);
	fill(C.begin(), C.end(), 0);

	//Select new sites at which to trigger storage
	vector < vector < int > > candidates = vector < vector < int > > (sites_pbwt_grouping.back() + 1);
	for (int l = 0 ; l < n_total_variants ; l++) if (sites_pbwt_evaluation[l]) candidates[sites_pbwt_grouping[l]].push_back(l);
	sites_pbwt_selection = vector < bool > (n_site , false);
	for (int g = 0 ; g < candidates.size() ; g++) {
		if (candidates[g].size() > 0) {
			sites_pbwt_selection[candidates[g][rng.getInt(candidates[g].size())]] = true;
		}
	}

	//Clean up previous selected states
	indexes_pbwt_neighbour = vector < vector < unsigned int > > (n_samples);

	//PBWT sweep
	for (int l = 0 ; l < n_total_variants ; l ++) {
		bool eval = sites_pbwt_evaluation[l];
		bool selc = sites_pbwt_selection[l];

		if (eval) {
			int u = 0, v = 0, p = l, q = l;
			for (int h = 0 ; h < n_haplotypes ; h ++) {
				int alookup = A[h], dlookup = C[h];
				if (dlookup > p) p = dlookup;
				if (dlookup > q) q = dlookup;
				if (!H_opt_var.get(l, alookup)) {
					A[u] = alookup;
					C[u] = p;
					p = 0;
					u++;
				} else {
					B[v] = alookup;
					D[v] = q;
					q = 0;
					v++;
				}
			}
			std::copy(B.begin(), B.begin()+v, A.begin()+u);
			std::copy(D.begin(), D.begin()+v, C.begin()+u);
			if (selc) store(l, A, C, M);
		}

		vrb.progress("  * PBWT section", l * 1.0 / n_total_variants);
	}
	vrb.bullet("PBWT selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void conditioning_set::store(int l, vector < int > & A, vector < int > & C, vector < int > & M) {
	for (int h = 0 ; h < n_haplotypes ; h ++) {
		int chap = A[h], cind = chap / 2, add_guess0 = 0, add_guess1 = 0, offset0 = 1, offset1 = 1, hap_guess0 = -1, hap_guess1 = -1, div_guess0 = -1, div_guess1 = -1;
		for (int n_added = 0 ; n_added < depth ; ) {
			if ((h-offset0)>=0) {
				hap_guess0 = A[h-offset0];
				div_guess0 = max(C[h-offset0+1], div_guess0);
				add_guess0 = 1;
			} else { add_guess0 = 0; div_guess0 = l+1; }
			if ((h+offset1)<n_hap) {
				hap_guess1 = A[h+offset1];
				div_guess1 = max(C[h+offset1], div_guess1);
				add_guess1 = 1;
			} else { add_guess1 = 0; div_guess1 = l+1; }
			if (add_guess0 && add_guess1) {
				if (div_guess0 < div_guess1) {
					if (hap_guess0 != M[chap * depth + n_added]) {
						indexes_pbwt_neighbour[cind].push_back(hap_guess0);
						M[chap * depth + n_added] = hap_guess0;
					}
					offset0++; n_added++;
				} else {
					if (hap_guess1 != M[chap * depth + n_added]) {
						indexes_pbwt_neighbour[cind].push_back(hap_guess1);
						M[chap * depth + n_added] = hap_guess1;
					}
					offset1++; n_added++;
				}
			} else if (add_guess0) {
				if (hap_guess0 != M[chap * depth + n_added]) {
					indexes_pbwt_neighbour[cind].push_back(hap_guess0);
					M[chap * depth + n_added] = hap_guess0;
				}
				offset0++; n_added++;
			} else if (add_guess1) {
				if (hap_guess1 != M[chap * depth + n_added]) {
					indexes_pbwt_neighbour[cind].push_back(hap_guess1);
					M[chap * depth + n_added] = hap_guess1;
				}
				offset1++; n_added++;
			} else {
				offset0++;
				offset1++;
			}
		}
	}
}
