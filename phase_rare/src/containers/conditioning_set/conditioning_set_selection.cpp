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

void conditioning_set::select(variant_map & V, genotype_set & G) {
	tac.clock();

	npushes = 0;
	ncollisions = 0;

	vector < int > A = vector < int > (n_haplotypes, 0);
	vector < int > B = vector < int > (n_haplotypes, 0);
	vector < int > R = vector < int > (n_haplotypes, 0);
	vector < int > M = vector < int > (depth_common * n_haplotypes, -1);
	iota(A.begin(), A.end(), 0);
	random_shuffle(A.begin(), A.end());

	//Select new sites at which to trigger storage
	vector < vector < int > > candidates = vector < vector < int > > (sites_pbwt_grouping.back() + 1);
	for (int l = 0 ; l < n_scaffold_variants ; l++) if (sites_pbwt_evaluation[l]) candidates[sites_pbwt_grouping[l]].push_back(l);
	sites_pbwt_selection = vector < bool > (n_scaffold_variants , false);
	for (int g = 0 ; g < candidates.size() ; g++) {
		if (candidates[g].size() > 0) {
			sites_pbwt_selection[candidates[g][rng.getInt(candidates[g].size())]] = true;
		}
	}

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
				int u = 0, v = 0;
				for (int h = 0 ; h < n_haplotypes ; h ++) {
					if (!Hvar.get(vs, A[h])) A[u++] = A[h];
					else B[v++] = A[h];
				}
				std::copy(B.begin(), B.begin()+v, A.begin()+u);
				for (int h = 0 ; h < n_haplotypes ; h ++) R[A[h]] = h;
				if (selc) storeCommon(A, M);
			}
		} else if (vr >= 0 && G.GRvar_genotypes[vr].size() > 1) storeRare(R, G.GRvar_genotypes[vr]);
		vrb.progress("  * PBWT forward selection", vt * 1.0 / V.sizeFull());
	}
	sort(indexes_pbwt_neighbour_serialized.begin(), indexes_pbwt_neighbour_serialized.end());
	indexes_pbwt_neighbour_serialized.erase(unique(indexes_pbwt_neighbour_serialized.begin(), indexes_pbwt_neighbour_serialized.end()), indexes_pbwt_neighbour_serialized.end());
	vrb.bullet("PBWT forward selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	//PBWT backward sweep
	tac.clock();
	for (int vt = V.sizeFull()-1 ; vt >= 0 ; vt --) {
		int vc = V.vec_full[vt]->idx_common;
		int vr = V.vec_full[vt]->idx_rare;
		int vs = V.vec_full[vt]->idx_scaffold;

		if (vs >= 0) {
			bool eval = sites_pbwt_evaluation[vs];
			bool selc = sites_pbwt_selection[vs];
			if (eval) {
				int u = 0, v = 0;
				for (int h = 0 ; h < n_haplotypes ; h ++) {
					if (!Hvar.get(vs, A[h])) A[u++] = A[h];
					else B[v++] = A[h];
				}
				std::copy(B.begin(), B.begin()+v, A.begin()+u);
				for (int h = 0 ; h < n_haplotypes ; h ++) R[A[h]] = h;
				if (selc) storeCommon(A, M);
			}
		} else if (vr >= 0 && G.GRvar_genotypes[vr].size() > 1) storeRare(R, G.GRvar_genotypes[vr]);
		vrb.progress("  * PBWT backward selection", (V.sizeFull()-vt) * 1.0 / V.sizeFull());
	}
	sort(indexes_pbwt_neighbour_serialized.begin(), indexes_pbwt_neighbour_serialized.end());
	indexes_pbwt_neighbour_serialized.erase(unique(indexes_pbwt_neighbour_serialized.begin(), indexes_pbwt_neighbour_serialized.end()), indexes_pbwt_neighbour_serialized.end());
	vrb.bullet("PBWT backward selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	stats1D statK;
	for (long int h = 0, e = 0 ; h < n_haplotypes ; h ++) {
		vector < unsigned int > buffer;
		while (indexes_pbwt_neighbour_serialized[e].first == h) {
			buffer.push_back(indexes_pbwt_neighbour_serialized[e].second);
			e++;
		}

		//Minimal number of states is 50
		if (buffer.size() < 50) {
			while (buffer.size() < 50) {
				if (!checkIBD2(shuffledO[shuffledI], h)) buffer.push_back(rng.getInt(n_haplotypes));
				shuffledI = (shuffledI<(n_haplotypes-1))?(shuffledI+1):0;
			}
			sort(buffer.begin(), buffer.end());
			buffer.erase(unique(buffer.begin(), buffer.end()), buffer.end());
		}

		indexes_pbwt_neighbour[h].reserve(buffer.size());
		indexes_pbwt_neighbour[h] = buffer;
		assert(indexes_pbwt_neighbour[h].size());
		statK.push(indexes_pbwt_neighbour[h].size());
	}

	indexes_pbwt_neighbour_serialized.clear();
	indexes_pbwt_neighbour_serialized.shrink_to_fit();
	vrb.bullet2("#states="+ stb.str(statK.mean(), 2) + "+/-" + stb.str(statK.sd(), 2));
	vrb.bullet2("#collisions = "+ stb.str(ncollisions) + " / #pushes = "+ stb.str(npushes) + " / rate = " + stb.str(npushes * 100.0 / (npushes + ncollisions), 2) + "%");
}

void conditioning_set::storeRare(vector < int > & R, vector < sparse_genotype > & G) {
	vector < pair < int, int > > N;
	for (int g = 0 ; g < G.size() ; g ++) {
		if (!G[g].mis) {
			unsigned int hap0 = 2*G[g].idx+0;
			unsigned int hap1 = 2*G[g].idx+1;
			N.push_back( pair < int, int > (R[hap0], hap0));
			N.push_back( pair < int, int > (R[hap1], hap1));
		}
	}
	sort(N.begin(), N.end());

	for (int h = 0 ; h < N.size() ; h ++) {
		int target_hap = N[h].second;
		int nstored = 0, offset = 1, done = 0;
		while (nstored < (2*depth_rare) && !done) {
			done = 1;
			if (h-offset >= 0) {
				//if (N[h-offset].second/2 != target_hap/2) {
				if (!checkIBD2(N[h-offset].second, target_hap)) {
					indexes_pbwt_neighbour_serialized.push_back(pair < unsigned int, unsigned int > (target_hap, N[h-offset].second));
					nstored ++;
				}
				done = 0;
			}
			if (h+offset < N.size()) {
				//if (N[h+offset].second/2 != target_hap/2) {
				if (!checkIBD2(N[h+offset].second, target_hap)) {
					indexes_pbwt_neighbour_serialized.push_back(pair < unsigned int, unsigned int > (target_hap, N[h+offset].second));
					nstored ++;
				}
				done = 0;
			}
			offset ++;
		}
	}
}

void conditioning_set::storeCommon(vector < int > & A, vector < int > & M) {
	for (int h = 0 ; h < n_haplotypes ; h ++) {
		int chap = A[h], add_guess0 = 0, add_guess1 = 0, offset0 = 1, offset1 = 1, hap_guess0 = -1, hap_guess1 = -1;
		for (int n_added = 0 ; n_added < depth_common ; ) {
			if ((h-offset0)>=0) {
				hap_guess0 = A[h-offset0];
				//add_guess0 = (hap_guess0/2 != chap/2);
				add_guess0 = !checkIBD2(hap_guess0/2, chap/2);
			} else add_guess0 = 0;
			if ((h+offset1)<n_haplotypes) {
				hap_guess1 = A[h+offset1];
				//add_guess1 = (hap_guess1/2 != chap/2);
				add_guess1 = !checkIBD2(hap_guess1/2, chap/2);
			} else add_guess1 = 0;
			if (add_guess0 && add_guess1) {
				if (hap_guess0 != M[chap * depth_common + n_added]) {
					indexes_pbwt_neighbour_serialized.push_back(pair < unsigned int, unsigned int > (chap, hap_guess0));
					M[chap * depth_common + n_added] = hap_guess0;
					npushes++;
				} else ncollisions++;
				offset0++; n_added++;
				if (hap_guess1 != M[chap * depth_common + n_added]) {
					indexes_pbwt_neighbour_serialized.push_back(pair < unsigned int, unsigned int > (chap, hap_guess1));
					M[chap * depth_common + n_added] = hap_guess1;
					npushes++;
				} else ncollisions++;
				offset1++; n_added++;
			} else if (add_guess0) {
				if (hap_guess0 != M[chap * depth_common + n_added]) {
					indexes_pbwt_neighbour_serialized.push_back(pair < unsigned int, unsigned int > (chap, hap_guess0));
					M[chap * depth_common + n_added] = hap_guess0;
					npushes++;
				} else ncollisions++;
				offset0++; n_added++;
			} else if (add_guess1) {
				if (hap_guess1 != M[chap * depth_common + n_added]) {
					indexes_pbwt_neighbour_serialized.push_back(pair < unsigned int, unsigned int > (chap, hap_guess1));
					M[chap * depth_common + n_added] = hap_guess1;
					npushes++;
				} else ncollisions++;
				offset1++; n_added++;
			} else {
				offset0++;
				offset1++;
			}
		}
	}
}

