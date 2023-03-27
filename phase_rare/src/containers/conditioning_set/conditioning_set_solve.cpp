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

void conditioning_set::solve(variant_map & V, genotype_set & G) {
	tac.clock();

	//
	vector < int > A = vector < int > (n_haplotypes, 0);
	vector < int > B = vector < int > (n_haplotypes, 0);
	vector < int > C = vector < int > (n_haplotypes, 0);
	vector < int > D = vector < int > (n_haplotypes, 0);
	vector < int > R = vector < int > (n_haplotypes, 0);
	iota(A.begin(), A.end(), 0);
	shuffle(A.begin(), A.end(), rng.getEngine());

	//Get cM positions of the scaffold sites
	vector < float > vs_cm = vector < float > (n_scaffold_variants, 0.0);
	for (int l = 0 ; l < n_scaffold_variants ; ++l) vs_cm[l] = V.vec_scaffold[l]->cm;

	//PBWT forward sweep
	tac.clock();
	for (int vt = 0 ; vt < V.sizeFull() ; vt ++) {
		int vc = V.vec_full[vt]->idx_common;
		int vr = V.vec_full[vt]->idx_rare;
		int vs = V.vec_full[vt]->idx_scaffold;

		if (vs >= 0) {
			bool eval = sites_pbwt_evaluation[vs];
			bool selc = sites_pbwt_selection[vs];
			if (eval) {
				int u = 0, v = 0, p = vs, q = vs;
				for (int h = 0 ; h < n_haplotypes ; h ++) {
					int alookup = A[h], dlookup = C[h];
					if (dlookup > p) p = dlookup;
					if (dlookup > q) q = dlookup;
					if (!Hvar.get(vs, alookup)) {
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
				for (int h = 0 ; h < n_haplotypes ; h ++) R[A[h]] = h;
			}
		} else if (vr >= 0) solveRareForward(A, C, R, G, vr, V.vec_rare[vr]->cm, vs_cm);
		vrb.progress("  * PBWT forward pass", vt * 1.0 / V.sizeFull());
	}
	vrb.bullet("PBWT forward pass (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	//PBWT backward sweep
	tac.clock();
	fill(C.begin(), C.end(), V.sizeScaffold() - 1);
	for (int vt = V.sizeFull()-1 ; vt >= 0 ; vt --) {
		int vc = V.vec_full[vt]->idx_common;
		int vr = V.vec_full[vt]->idx_rare;
		int vs = V.vec_full[vt]->idx_scaffold;

		if (vs >= 0) {
			bool eval = sites_pbwt_evaluation[vs];
			bool selc = sites_pbwt_selection[vs];
			if (eval) {
				int u = 0, v = 0, p = vs, q = vs;
				for (int h = 0 ; h < n_haplotypes ; h ++) {
					int alookup = A[h], dlookup = C[h];
					if (dlookup < p) p = dlookup;
					if (dlookup < q) q = dlookup;
					if (!Hvar.get(vs, alookup)) {
						A[u] = alookup;
						C[u] = p;
						p = V.sizeScaffold() - 1;
						u++;
					} else {
						B[v] = alookup;
						D[v] = q;
						q = V.sizeScaffold() - 1;
						v++;
					}
				}
				std::copy(B.begin(), B.begin()+v, A.begin()+u);
				std::copy(D.begin(), D.begin()+v, C.begin()+u);
				for (int h = 0 ; h < n_haplotypes ; h ++) R[A[h]] = h;
			}
		} else if (vr >= 0) solveRareBackward(A, C, R, G, vr, V.vec_rare[vr]->cm, vs_cm);
		vrb.progress("  * PBWT backward pass", (V.sizeFull() - vt) * 1.0 / V.sizeFull());
	}
	vrb.bullet("PBWT backward pass (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	assert(!CF.size());
}

void conditioning_set::solveRareForward(vector < int > & A, vector < int > & D, vector < int > & R, genotype_set & G, unsigned int vr, float vr_cm, vector < float > & vs_cm) {
	vector < int > S, C = vector < int > (n_haplotypes, G.major_alleles[vr]?1:-1);
	for (int g = 0 ; g < G.GRvar_genotypes[vr].size() ; g ++) {
		if (G.GRvar_genotypes[vr][g].pha) {
			C[2*G.GRvar_genotypes[vr][g].idx+0] = G.GRvar_genotypes[vr][g].al0?1:-1;
			C[2*G.GRvar_genotypes[vr][g].idx+1] = G.GRvar_genotypes[vr][g].al1?1:-1;
		} else {
			C[2*G.GRvar_genotypes[vr][g].idx+0] = 0;
			C[2*G.GRvar_genotypes[vr][g].idx+1] = 0;
			S.push_back(g);
		}
	}

	//PHASING FIRST PASS
	float thresh = 2.5, v, v0, v1;
	while (S.size() && thresh > 1.0) {
		unsigned int sizeS = S.size();
		for (vector < int > :: iterator s = S.begin() ; s != S.end() ; ) {
			int h0 = G.GRvar_genotypes[vr][*s].idx*2+0;
			int h1 = G.GRvar_genotypes[vr][*s].idx*2+1;

			if (R[h0]>0) v0 = C[A[R[h0]-1]];
			if (R[h0]<(n_haplotypes-1)) v0 += C[A[R[h0]+1]];
			if (R[h1]>0) v1 = C[A[R[h1]-1]];
			if (R[h1]<(n_haplotypes-1)) v1 += C[A[R[h1]+1]];
			v = v0 - v1;

			if (v > thresh) {
				C[h0] = 1.0; C[h1] = -1.0;
				G.GRvar_genotypes[vr][*s].phase(2);
				s = S.erase(s);
			} else if (v < -thresh) {
				C[h0] = -1.0; C[h1] = 1.0;
				G.GRvar_genotypes[vr][*s].phase(1);
				s = S.erase(s);
			} else s++;
		}
		if (S.size() == sizeS) thresh -= 1.0 ;
	}

	//PHASING SECOND PASS
	for (vector < int > :: iterator s = S.begin() ; s != S.end() ; s++) {
		int h0 = G.GRvar_genotypes[vr][*s].idx*2+0;
		int h1 = G.GRvar_genotypes[vr][*s].idx*2+1;

		v0 = v1 = 0;
		if (R[h0]>0) v0 += C[A[R[h0]-1]] * abs(vr_cm - vs_cm[D[R[h0]]]);
		if (R[h0]<(n_haplotypes-1)) v0 += C[A[R[h0]+1]] * abs(vr_cm - vs_cm[D[R[h0]+1]]);
		if (R[h1]>0) v1 -= C[A[R[h1]-1]] * abs(vr_cm - vs_cm[D[R[h1]]]);
		if (R[h1]<(n_haplotypes-1)) v1 -= C[A[R[h1]+1]] * abs(vr_cm - vs_cm[D[R[h1]+1]]);
		v = v0+v1;

		if (v > 0) {
			C[h0] = 1; C[h1] = -1;
			CF.emplace_back(2, v);
		} else {
			C[h0] = -1; C[h1] = 1;
			CF.emplace_back(1, v);
		}
	}
}


void conditioning_set::solveRareBackward(vector < int > & A, vector < int > & D, vector < int > & R, genotype_set & G, unsigned int vr, float vr_cm, vector < float > & vs_cm) {
	vector < int > S, C = vector < int > (n_haplotypes, G.major_alleles[vr]?1:-1);
	for (int g = G.GRvar_genotypes[vr].size()-1 ; g >= 0 ; g --) {
		if (G.GRvar_genotypes[vr][g].pha) {
			C[2*G.GRvar_genotypes[vr][g].idx+0] = G.GRvar_genotypes[vr][g].al0?1:-1;
			C[2*G.GRvar_genotypes[vr][g].idx+1] = G.GRvar_genotypes[vr][g].al1?1:-1;
		} else {
			C[2*G.GRvar_genotypes[vr][g].idx+0] = 0;
			C[2*G.GRvar_genotypes[vr][g].idx+1] = 0;
			S.push_back(g);
		}
	}

	//PHASING FIRST PASS
	vector < int > Stmp = S;
	float thresh = 2.5, v, v0, v1;
	while (Stmp.size() && thresh > 1.0) {
		unsigned int sizeS = Stmp.size();
		for (vector < int > :: iterator s = Stmp.begin() ; s != Stmp.end() ; ) {
			int h0 = G.GRvar_genotypes[vr][*s].idx*2+0;
			int h1 = G.GRvar_genotypes[vr][*s].idx*2+1;

			if (R[h0]>0) v0 = C[A[R[h0]-1]];
			if (R[h0]<(n_haplotypes-1)) v0 += C[A[R[h0]+1]];
			if (R[h1]>0) v1 = C[A[R[h1]-1]];
			if (R[h1]<(n_haplotypes-1)) v1 += C[A[R[h1]+1]];
			v = v0 - v1;

			if (v > thresh) {
				C[h0] = 1.0; C[h1] = -1.0;
				G.GRvar_genotypes[vr][*s].phase(2);
				s = Stmp.erase(s);
			} else if (v < -thresh) {
				C[h0] = -1.0; C[h1] = 1.0;
				G.GRvar_genotypes[vr][*s].phase(1);
				s = Stmp.erase(s);
			} else s++;
		}
		if (Stmp.size() == sizeS) thresh -= 1.0 ;
	}

	//PHASING SECOND PASS
	cflip ctmp;
	for (vector < int > :: iterator s = S.begin() ; s != S.end() ; s++) {
		int h0 = G.GRvar_genotypes[vr][*s].idx*2+0;
		int h1 = G.GRvar_genotypes[vr][*s].idx*2+1;

		if (find(Stmp.begin(), Stmp.end(), *s)!=Stmp.end()) {
			v0 = v1 = 0;
			if (R[h0]>0) v0 += C[A[R[h0]-1]] * abs(vr_cm - vs_cm[D[R[h0]]]);
			if (R[h0]<(n_haplotypes-1)) v0 += C[A[R[h0]+1]] * abs(vr_cm - vs_cm[D[R[h0]+1]]);
			if (R[h1]>0) v1 -= C[A[R[h1]-1]] * abs(vr_cm - vs_cm[D[R[h1]]]);
			if (R[h1]<(n_haplotypes-1)) v1 -= C[A[R[h1]+1]] * abs(vr_cm - vs_cm[D[R[h1]+1]]);
			v = v0+v1;

			if (v > 0) {
				C[h0] = 1; C[h1] = -1;
				ctmp.set(2, v);
			} else {
				C[h0] = -1; C[h1] = 1;
				ctmp.set(1, v);
			}

			G.GRvar_genotypes[vr][*s].pha = 1;
			if (!ctmp.betterThan(CF.back())) ctmp = CF.back();

			if (ctmp.pgenotype == 1) {
				G.GRvar_genotypes[vr][*s].al0 = 0;
				G.GRvar_genotypes[vr][*s].al1 = 1;
			} else {
				G.GRvar_genotypes[vr][*s].al0 = 1;
				G.GRvar_genotypes[vr][*s].al1 = 0;
			}
		}
		CF.pop_back();
	}
}
