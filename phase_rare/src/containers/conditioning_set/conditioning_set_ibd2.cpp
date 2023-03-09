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

#include <containers/conditioning_set/conditioning_set_header.h>

using namespace std;

void conditioning_set::scanIBD2(variant_map & V) {
	tac.clock();

	//
	int M = 3, n_ind = n_haplotypes / 2;
	vector < int > U = vector < int > (M, 0);
	vector < int > P = vector < int > (M, 0);
	vector < int > G = vector < int > (n_ind, 0);
	vector < vector < int > > A = vector < vector < int > > (M, vector < int > (n_ind, 0));
	vector < vector < int > > D = vector < vector < int > > (M, vector < int > (n_ind, 0));

	//
	for (int l = 0 ; l < n_scaffold_variants ; l ++) {
		fill(U.begin(), U.end(), 0);
		fill(P.begin(), P.end(), l);
		for (int i = 0 ; i < n_ind ; i ++) {
			int alookup = l?A[0][i]:i;
			int dlookup = l?D[0][i]:0;
			for (int g = 0 ; g < M ; g++) if (dlookup > P[g]) P[g] = dlookup;
			G[i] = Hvar.get(l, 2*alookup+0) + Hvar.get(l, 2*alookup+1);
			A[G[i]][U[G[i]]] = alookup;
			D[G[i]][U[G[i]]] = P[G[i]];
			P[G[i]] = 0;
			U[G[i]]++;
		}
		for (int g = 1, offset = U[0] ; g < M ; g++) {
			copy(A[g].begin(), A[g].begin()+U[g], A[0].begin() + offset);
			copy(D[g].begin(), D[g].begin()+U[g], D[0].begin() + offset);
			offset += U[g];
		}
		for (int i = 1 ; i < n_ind ; i ++) {
			int ind0 = A[0][i];
			int ng0 = (l<(n_scaffold_variants-1))?(Hvar.get(l+1, 2*ind0+0)+Hvar.get(l+1, 2*ind0+1)):-1;
			for (int ip = i-1, div = -1 ; ip >= 0 ; ip --) {
				if (G[ip] != G[i]) break;
				div = max(div, D[0][ip+1]);
				double lengthMatchCM = V.vec_scaffold[l]->cm - V.vec_scaffold[div]->cm;
				double lengthMatchBP = V.vec_scaffold[l]->bp - V.vec_scaffold[div]->bp;
				double lengthMatchCT = l - div + 1;
				if ((lengthMatchCT == n_scaffold_variants) || (lengthMatchCM >= 2.5f && lengthMatchBP >= 1e6 && lengthMatchCT >= 100)) {
					int ind1 = A[0][ip];
					int ng1 = (l<(n_scaffold_variants-1))?(Hvar.get(l+1, 2*ind1+0)+Hvar.get(l+1, 2*ind1+1)):-1;
					if (ng0 < 0 || ng0 != ng1) {
						//int disc = 0;
						//for (int hh = div ; hh <= l ; hh ++) disc += ((Hvar.get(hh, 2*ind1+0)+Hvar.get(hh, 2*ind1+1)) - (Hvar.get(hh, 2*ind0+0)+Hvar.get(hh, 2*ind0+1)));
						//cout << "IBD2: " << ind0 << " " << ind1 << " " << lengthMatchCM << " " << lengthMatchBP << " " << lengthMatchCT << " " << disc << endl;
						IBD2[min(ind0, ind1)].push_back(max(ind0, ind1));
					}
				} else break;
			}
		}
		vrb.progress("  * IBD2 constraints ", (l+1)*1.0/n_scaffold_variants);
	}

	unsigned long npairstot = 0, npairsind = 0;
	for (int i = 0 ; i < n_ind ; i ++) {
		sort(IBD2[i].begin(), IBD2[i].end());
		IBD2[i].erase(unique(IBD2[i].begin(), IBD2[i].end()), IBD2[i].end());
		npairstot += IBD2[i].size();
		npairsind += (IBD2[i].size()>0);
	}
	vrb.bullet("IBD2 constraints [#inds=" + stb.str(npairsind) + " / #pairs=" + stb.str(npairstot) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
