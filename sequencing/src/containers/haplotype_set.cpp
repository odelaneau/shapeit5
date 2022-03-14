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
	n_common_variants = 0.0;
	n_rare_variants = 0.0;
	n_total_variants = 0.0;
	n_total_haplotypes = 0.0;
	n_target_haplotypes = 0.0;
	n_reference_haplotypes = 0.0;
	n_total_samples = 0.0;
	n_target_samples = 0.0;
	n_reference_samples = 0.0;
	pbwt_modulo = 0.0;
	pbwt_depth = 0;
	pbwt_mac = 0;
	pbwt_mdr = 0.0;
	pbwt_nstored = 0;
	nthreads = 0;
	pbwt_evaluated.clear();
	pbwt_stored.clear();
	pbwt_parray.clear();
	pbwt_darray.clear();
	pbwt_neighbours.clear();
}

void haplotype_set::allocate(unsigned int _n_target_samples, unsigned int _n_reference_samples, unsigned int _n_common_variants, unsigned int _n_rare_variants, vector < bool > & _flag_common) {
	tac.clock();
	flag_common = _flag_common;
	n_target_samples = _n_target_samples;
	n_reference_samples = _n_reference_samples;
	n_total_samples = n_reference_samples + n_target_samples;
	n_target_haplotypes = 2 * n_target_samples;
	n_reference_haplotypes =  2 * n_reference_samples;
	n_total_haplotypes = 2 * n_total_samples;
	n_common_variants = _n_common_variants;
	n_rare_variants = _n_rare_variants;
	n_total_variants = n_common_variants + n_rare_variants;

	Hvar.allocate(n_common_variants, n_total_haplotypes);
	Hhap.allocate(n_total_haplotypes, n_common_variants);

	if (n_rare_variants>0) {
		Svar = vector < vector < unsigned int > > (n_total_variants, vector < unsigned int > ());
		Shap = vector < vector < unsigned int > > (n_total_haplotypes, vector < unsigned int > ());
	}

	vrb.bullet("HAP allocation [n_variants=" + stb.str(n_common_variants) + ", " + stb.str(n_rare_variants) + " / n_haplotypes = " + stb.str(n_target_samples) + ", " + stb.str(n_reference_samples) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::parametrizePBWT(int _pbwt_depth, double _pbwt_modulo, int _pbwt_mac, double _pbwt_mdr, int _nthreads) {
	pbwt_modulo = _pbwt_modulo;
	pbwt_depth = _pbwt_depth;
	pbwt_mac = _pbwt_mac;
	pbwt_mdr = _pbwt_mdr;
	nthreads = _nthreads;
}

void haplotype_set::initializePBWTmapping(variant_map & V) {
	tac.clock();
	for (int lc = 0 ; lc < n_common_variants ; lc ++) {
		if (V.vec_common[lc]->getMAC() >= pbwt_mac && V.vec_common[lc]->getMDR() <= pbwt_mdr) {
			pbwt_evaluated.push_back(lc);
			pbwt_cm.push_back(V.vec_common[lc]->cm);
			pbwt_grp.push_back((int)round(V.vec_common[lc]->cm / pbwt_modulo));
		}
	}
	for (int l = 0, src = -1, tar = -1 ; l < pbwt_grp.size() ; l ++) {
		if (src == pbwt_grp[l]) pbwt_grp[l] = tar;
		else {
			src = pbwt_grp[l];
			pbwt_grp[l] = ++tar;
		}
	}
	pbwt_nstored = pbwt_grp.back() + 1;
	updatePBWTmapping();
	vrb.bullet("PBWT indexing [l=" + stb.str(pbwt_nstored) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::updatePBWTmapping() {
	pbwt_stored = vector < int > (pbwt_grp.size() , -1);
	for (int idx = 0, loffset = 0, liter = 0 ; idx <= pbwt_grp.back() ; idx ++) {
		for (liter = 0; (loffset+liter) < pbwt_grp.size() && (pbwt_grp[loffset+liter]==idx) ; liter++);
		pbwt_stored[loffset + rng.getInt(liter)] = idx;
		loffset += liter;
	}
}

void haplotype_set::allocatePBWTarrays() {
	assert(pbwt_evaluated.size() > 0);
	//pbwt_neighbours = vector < int > ((pbwt_depth+1) * pbwt_nstored * n_target_samples * 2UL, 0);	//The depth+1 is there for transposing the matrix
	pbwt_neighbours = vector < vector < unsigned int > > (n_target_samples);	//The depth+1 is there for transposing the matrix
	pbwt_parray = vector < int > (n_total_haplotypes, 0);
	pbwt_darray = vector < int > (n_total_haplotypes, 0);
}

void haplotype_set::updateHaplotypes(genotype_set & G, bool first_time) {
	tac.clock();
	for (unsigned int i = 0 ; i < G.n_target_samples ; i ++) {

		//COMMONS
		for (unsigned int vc = 0 ; vc < n_common_variants ; vc ++) {
			if (first_time || (VAR_GET_HET(MOD2(vc), G.vecG[i]->Variants[DIV2(vc)])) || (VAR_GET_MIS(MOD2(vc), G.vecG[i]->Variants[DIV2(vc)]))) {
				bool a0 = VAR_GET_HAP0(MOD2(vc), G.vecG[i]->Variants[DIV2(vc)]);
				bool a1 = VAR_GET_HAP1(MOD2(vc), G.vecG[i]->Variants[DIV2(vc)]);
				Hhap.set(2*i+0, vc, a0);
				Hhap.set(2*i+1, vc, a1);
			}
		}

		//RARE
		if (n_rare_variants > 0) {
			Shap[2*i+0].clear();
			Shap[2*i+1].clear();
			for (unsigned int vr = 0 ; vr < G.vecG[i]->RareIndexes.size() ; vr ++) {
				bool a0 = VAR_GET_HAP0(MOD2(vr), G.vecG[i]->RareVariants[DIV2(vr)]);
				bool a1 = VAR_GET_HAP1(MOD2(vr), G.vecG[i]->RareVariants[DIV2(vr)]);
				if (a0) Shap[2*i+0].push_back(G.vecG[i]->RareIndexes[vr]);
				if (a1) Shap[2*i+1].push_back(G.vecG[i]->RareIndexes[vr]);
				//cout << a0 << " " << a1 << endl;
			}
		}

	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_H2V(bool full) {
	tac.clock();

	//COMMON
	if (!full) Hhap.transpose(Hvar, 2*n_target_samples, n_common_variants);
	else Hhap.transpose(Hvar, n_total_haplotypes, n_common_variants);

	//RARE [To be optimize for non full transpose, i.e. split sparse for ref and target]
	if (n_rare_variants > 0) {
		for (unsigned int vt = 0 ; vt < n_total_variants ; vt ++) Svar[vt].clear();
		for (unsigned int h = 0 ; h < n_total_haplotypes ; h ++) for (unsigned int r = 0 ; r < Shap[h].size() ; r ++) Svar[Shap[h][r]].push_back(h);
	}

	vrb.bullet("H2V transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_V2H(bool full) {
	tac.clock();

	//COMMON
	if (!full) Hvar.transpose(Hhap, n_common_variants, 2 * n_target_samples);
	else Hvar.transpose(Hhap, n_common_variants, n_total_haplotypes);

	//RARE
	if (n_rare_variants > 0) {
		for (unsigned int h = 0 ; h < n_total_haplotypes ; h ++) Shap[h].clear();
		for (unsigned int vt = 0 ; vt < n_total_variants ; vt ++) if (!flag_common [vt]) for (unsigned int r = 0 ; r < Svar[vt].size() ; r ++) Shap[Svar[vt][r]].push_back(vt);
	}

	vrb.bullet("V2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
/*
void haplotype_set::transposePBWTarrays() {
	tac.clock();
	int block = 32;
	unsigned long addr_tar, addr_src;
	unsigned long addr_offset = pbwt_nstored * n_target_samples * 2UL;
	for (int d = 0; d < pbwt_depth ; d ++) {
		for (int s = 0; s < pbwt_nstored ; s += block) {
			for(int h = 0; h < n_target_samples * 2; ++h) {
				for(int b = 0; b < block && s + b < pbwt_nstored; ++b) {
					addr_tar = pbwt_depth * addr_offset + h*pbwt_nstored + s + b;
					addr_src = d * addr_offset + (s + b)*n_target_samples*2UL + h;
					pbwt_neighbours[addr_tar] = pbwt_neighbours[addr_src];
				}
			}
		}
		std::copy(pbwt_neighbours.begin() + pbwt_depth * addr_offset , pbwt_neighbours.end(), pbwt_neighbours.begin() + d * addr_offset );
	}
	vrb.bullet("C2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
*/
void haplotype_set::selectPBWTarrays() {
	tac.clock();
	if (bannedPairs.size() == 0) {
		bannedPairs = vector < vector < unsigned int > > (n_target_samples);
	}
	vector < int > B = vector < int > (n_total_haplotypes, 0);
	vector < int > D = vector < int > (n_total_haplotypes, 0);
	vector < int > M = vector < int > (n_target_haplotypes * pbwt_depth, -1);
	for (int i = 0 ; i < n_target_samples ; i ++) pbwt_neighbours[i].clear();
	for (int h = 0 ; h < n_total_haplotypes ; h ++) pbwt_parray[h] = h;
	fill(pbwt_darray.begin(), pbwt_darray.end(), 0);
	for (int l = 0 ; l < pbwt_evaluated.size() ; l ++) {
		int u = 0, v = 0, p = l, q = l;

		//PBWT PASS
		for (int h = 0 ; h < n_total_haplotypes ; h ++) {
			int alookup = pbwt_parray[h];
			int dlookup = pbwt_darray[h];
			if (dlookup > p) p = dlookup;
			if (dlookup > q) q = dlookup;
			if (!Hvar.get(pbwt_evaluated[l], alookup)) {
				pbwt_parray[u] = alookup;
				pbwt_darray[u] = p;
				p = 0;
				u++;
			} else {
				B[v] = alookup;
				D[v] = q;
				q = 0;
				v++;
			}
		}
		std::copy(B.begin(), B.begin()+v, pbwt_parray.begin()+u);
		std::copy(D.begin(), D.begin()+v, pbwt_darray.begin()+u);

		//PBWT STORAGE
		if (pbwt_stored[l] >= 0) {
			unsigned long addr_offset = pbwt_nstored * n_target_samples * 2UL;
			for (int h = 0 ; h < n_total_haplotypes ; h ++) {
				int chap = pbwt_parray[h];
				int cind = chap / 2;
				if (cind < n_target_samples) {
					int add_guess0 = 0, add_guess1 = 0, offset0 = 1, offset1 = 1, hap_guess0 = -1, hap_guess1 = -1, div_guess0 = -1, div_guess1 = -1;
					unsigned long tar_idx = pbwt_stored[l] * 2UL * n_target_samples + chap;
					for (int n_added = 0 ; n_added < pbwt_depth ; ) {
						if ((h-offset0)>=0) {
							hap_guess0 = pbwt_parray[h-offset0];
							div_guess0 = max(pbwt_darray[h-offset0+1], div_guess0);
							add_guess0 = checkIBD2matching(chap, hap_guess0);
						} else { add_guess0 = 0; div_guess0 = l+1; }
						if ((h+offset1)<n_total_haplotypes) {
							hap_guess1 = pbwt_parray[h+offset1];
							div_guess1 = max(pbwt_darray[h+offset1], div_guess1);
							add_guess1 = checkIBD2matching(chap, hap_guess1);
						} else { add_guess1 = 0; div_guess1 = l+1; }
						if (add_guess0 && add_guess1) {
							if (div_guess0 < div_guess1) {
								//pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess0;
								if (hap_guess0 != M[chap * pbwt_depth + n_added]) {
									M[chap * pbwt_depth + n_added] = hap_guess0;
									pbwt_neighbours[cind].push_back(hap_guess0);
								}
								offset0++; n_added++;
							} else {
								//pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess1;
								if (hap_guess1 != M[chap * pbwt_depth + n_added]) {
									M[chap * pbwt_depth + n_added] = hap_guess1;
									pbwt_neighbours[cind].push_back(hap_guess1);
								}
								offset1++; n_added++;
							}
						} else if (add_guess0) {
							//pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess0;
							if (hap_guess0 != M[chap * pbwt_depth + n_added]) {
								M[chap * pbwt_depth + n_added] = hap_guess0;
								pbwt_neighbours[cind].push_back(hap_guess0);
							}
							offset0++; n_added++;
						} else if (add_guess1) {
							//pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess1;
							if (hap_guess1 != M[chap * pbwt_depth + n_added]) {
								M[chap * pbwt_depth + n_added] = hap_guess1;
								pbwt_neighbours[cind].push_back(hap_guess1);
							}
							offset1++; n_added++;
						} else {
							offset0++;
							offset1++;
						}
					}
				}
			}
		}
		vrb.progress("  * PBWT selection", (l+1)*1.0/pbwt_evaluated.size());
	}
	vrb.bullet("PBWT selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::mergeIBD2constraints() {
	unsigned int n_ibd2_inds = 0;
	unsigned int n_ibd2_pairs = 0;
	for (int i = 0 ; i < n_target_samples ; i ++) {
		sort(bannedPairs[i].begin(), bannedPairs[i].end());
		bannedPairs[i].erase(unique(bannedPairs[i].begin(), bannedPairs[i].end()), bannedPairs[i].end());
		if (bannedPairs[i].size() > 0) n_ibd2_inds++;
		n_ibd2_pairs+=bannedPairs[i].size();
	}
	vrb.bullet("IBD2 constraints [#inds=" + stb.str(n_ibd2_inds) + " / #constraints=" + stb.str(n_ibd2_pairs) + "]");
}
