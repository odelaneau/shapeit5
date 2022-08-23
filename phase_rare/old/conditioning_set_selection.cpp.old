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

void conditioning_set::select2(variant_map & V, genotype_set & G) {
	tac.clock();

	vector < int > A = vector < int > (n_haplotypes, 0);
	vector < int > B = vector < int > (n_haplotypes, 0);
	vector < int > R = vector < int > (n_haplotypes, 0);
	iota(A.begin(), A.end(), 0);
	random_shuffle(A.begin(), A.end());

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
			}
		} else if (vr >= 0 && G.GRvar_genotypes[vr].size() > 1) {
			storeAtRareUnphased(R, G.GRvar_genotypes[vr]);
		} else if (vc >= 0) {
			storeAtCommonUnphased(R, G, vc);
		}
		vrb.progress("  * PBWT forward selection", vt * 1.0 / V.sizeFull());
	}
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
			}
		} else if (vr >= 0 && G.GRvar_genotypes[vr].size() > 1) {
			storeAtRareUnphased(R, G.GRvar_genotypes[vr]);
		} else if (vc >= 0) {
			storeAtCommonUnphased(R, G, vc);
		}
		vrb.progress("  * PBWT backward selection", (V.sizeFull()-vt) * 1.0 / V.sizeFull());
	}
	vrb.bullet("PBWT backward selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	basic_stats statK;
	for (int h = 0 ; h < n_haplotypes ; h ++) {
		sort(indexes_pbwt_neighbour[h].begin(), indexes_pbwt_neighbour[h].end());
		indexes_pbwt_neighbour[h].erase(unique(indexes_pbwt_neighbour[h].begin(), indexes_pbwt_neighbour[h].end()), indexes_pbwt_neighbour[h].end());

		assert(indexes_pbwt_neighbour[h].size());
		statK.push(indexes_pbwt_neighbour[h].size());
	}

	//indexes_pbwt_neighbour_serialized.clear();
	indexes_pbwt_neighbour_serialized.shrink_to_fit();
	vrb.bullet2("#states="+ stb.str(statK.mean(), 2) + "+/-" + stb.str(statK.sd(), 2));
}

struct rinfo {
	unsigned int hidx;
	bool phased;

	bool operator < (const rinfo & s) const {
		return (hidx<s.hidx);
	}

};

void conditioning_set::storeAtRareUnphased(vector < int > & R, vector < rare_genotype > & G) {
	vector < pair < int, rinfo > > N;
	for (int g = 0 ; g < G.size() ; g ++) {
		if (!G[g].mis) {
			rinfo r0, r1;
			r0.hidx = 2*G[g].idx+0;
			r1.hidx = 2*G[g].idx+1;
			r0.phased = G[g].pha;
			r1.phased = G[g].pha;
			N.push_back( pair < int, rinfo > (R[r0.hidx], r0));
			N.push_back( pair < int, rinfo > (R[r1.hidx], r1));
		}
	}

	sort(N.begin(), N.end());

	for (int h = 0 ; h < N.size() ; h ++) {
		bool targ_phased = N[h].second.phased;
		int targ_hap = N[h].second.hidx;

		int from_idx = (h>0)?(h-1):h;
		int from_cnt = 0;
		while (from_idx>=0 && from_cnt < depth_rare) {
			bool cond_phased = N[from_idx].second.phased;
			int cond_hap = N[from_idx].second.hidx;
			if (cond_phased && targ_hap/2 != cond_hap / 2) {
				indexes_pbwt_neighbour[targ_hap].push_back(cond_hap);
				from_cnt++;
			}
			from_idx--;
		}

		int to_idx = (h<(N.size()-1))?(h+1):h;
		int to_cnt = 0;
		while (to_idx<N.size() && to_cnt < depth_rare) {
			bool cond_phased = N[to_idx].second.phased;
			int cond_hap = N[to_idx].second.hidx;
			if (cond_phased && targ_hap/2 != cond_hap / 2) {
				indexes_pbwt_neighbour[targ_hap].push_back(cond_hap);
				to_cnt++;
			}
			to_idx++;
		}
	}
}

void conditioning_set::scrollUntilPhasedIsFound(vector < int > & A, genotype_set & G, int vc, int thap, int tidx, bool a) {
	int add_guess0 = 0, add_guess1 = 0, offset0 = 1, offset1 = 1, hap_guess0 = -1, hap_guess1 = -1;
	for (int n_added = 0 ; n_added < depth_common ; ) {

		if ((tidx-offset0)>=0) {
			hap_guess0 = A[tidx-offset0];
			add_guess0 = (hap_guess0/2 != thap/2) && G.GCvar_phased.get(vc, hap_guess0/2) && (G.GCvar_alleles.get(vc, hap_guess0) == a);
		} else add_guess0 = 0;
		if ((tidx+offset1)<n_haplotypes) {
			hap_guess1 = A[tidx+offset1];
			add_guess1 = (hap_guess1/2 != thap/2) && G.GCvar_phased.get(vc, hap_guess1/2) && (G.GCvar_alleles.get(vc, hap_guess1) == a);
		} else add_guess1 = 0;

		if (add_guess0 && add_guess1) {
			indexes_pbwt_neighbour[thap].push_back(hap_guess0);
			indexes_pbwt_neighbour[thap].push_back(hap_guess1);
			offset0++;
			offset1++;
			n_added+=2;
		} else if (add_guess0) {
			indexes_pbwt_neighbour[thap].push_back(hap_guess0);
			offset0++;
			n_added++;
		} else if (add_guess1) {
			indexes_pbwt_neighbour[thap].push_back(hap_guess1);
			offset1++;
			n_added++;
		} else {
			offset0++;
			offset1++;
		}
	}
}


void conditioning_set::storeAtCommonUnphased(vector < int > & A, genotype_set & G, int vc) {
	for (int h = 0 ; h < n_haplotypes ; h ++) {
		int thap = A[h];
		int tind = thap/2;

		if (!G.GCvar_phased.get(vc, tind)) {
			scrollUntilPhasedIsFound(A, G, vc, thap, h, false);
			scrollUntilPhasedIsFound(A, G, vc, thap, h, true);
		}
	}
}
