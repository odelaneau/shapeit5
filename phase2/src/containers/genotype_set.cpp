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
	n_samples = 0.0;
	names.clear();
	GRvar_genotypes.clear();
	GRind_genotypes.clear();
}

void genotype_set::imputeMonomorphic() {
	unsigned int n_solved = 0;
	for (int vr = 0 ; vr < n_rare_variants ; vr ++) {
		bool mono = true;
		for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {
			if (!GRvar_genotypes[vr][r].mis) mono = false;
		}
		if (mono) {
			for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {
				GRvar_genotypes[vr][r].al0 = major_alleles[vr];
				GRvar_genotypes[vr][r].al1 = major_alleles[vr];
				GRvar_genotypes[vr][r].pha = 1;
				n_solved++;
			}
		}
	}
	vrb.bullet(stb.str(n_solved) + " missing genotypes imputed at monomorphic sites");
}

void genotype_set::randomizeSingleton() {
	unsigned int n_solved = 0;
	for (int vr = 0 ; vr < n_rare_variants ; vr ++) {
		if (GRvar_genotypes[vr].size() == 1 && GRvar_genotypes[vr][0].het) {
			if (rng.flipCoin()) {
				GRvar_genotypes[vr][0].al0 = 0;
				GRvar_genotypes[vr][0].al1 = 1;
			} else {
				GRvar_genotypes[vr][0].al0 = 1;
				GRvar_genotypes[vr][0].al1 = 0;
			}
			GRvar_genotypes[vr][0].pha = 1;
			n_solved++;
		}
	}
	vrb.bullet(stb.str(n_solved) + " singletons randomly phased");
}

void genotype_set::impute(unsigned int vr, unsigned int hap, vector < float > & alphaXbeta_prev, vector < float > & alphaXbeta_curr, vector < unsigned int > & H) {
	float p[2] = { 0.0f };

	int tidx = -1;
	for (int k = 0, r = 0; (k< H.size()) || (r<GRvar_genotypes[vr].size()) ;) {
		int tind = (r<GRvar_genotypes[vr].size())?GRvar_genotypes[vr][r].idx:-1;
		int tmis = (r<GRvar_genotypes[vr].size())?GRvar_genotypes[vr][r].mis:-1;
		int cind = (k<H.size())?H[k]/2:-1;

		if (tind == hap/2) {
			tidx = r;
			r++;
		} else if (tind < 0) {
			p[major_alleles[vr]] += alphaXbeta_prev[k] * 0.5f + alphaXbeta_curr[k] * 0.5f;
			k++;
		} else if (cind < 0) {
			r++;
		} else if (tind < cind) {
			r++;
		} else if (tind==cind) {
			if (!tmis) p[!major_alleles[vr]] += alphaXbeta_prev[k] * 0.5f + alphaXbeta_curr[k] * 0.5f;
			k++;
		} else {
			p[major_alleles[vr]] += alphaXbeta_prev[k] * 0.5f + alphaXbeta_curr[k] * 0.5f;
			k++;
		}
	}

	assert (tidx>=0);
	if (hap%2) GRvar_genotypes[vr][tidx].ph1 = p[1] / (p[0] + p[1]);
	else GRvar_genotypes[vr][tidx].ph0 = p[1] / (p[0] + p[1]);
}

void genotype_set::phase(unsigned long int & n_phased, unsigned long int & n_total, float threshold) {
	n_phased = n_total = 0;

	float ee = 0.9999, ed = 0.0001;
	vector < float > gprobs = vector < float > (4, 0.0f);
	for (int vr = 0 ; vr < GRvar_genotypes.size() ; vr ++) {
		for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {

			if (GRvar_genotypes[vr][r].het) {
				float p00 = 1.0f - GRvar_genotypes[vr][r].ph0;
				float p01 = GRvar_genotypes[vr][r].ph0;
				float p10 = 1.0f - GRvar_genotypes[vr][r].ph1;
				float p11 = GRvar_genotypes[vr][r].ph1;
				float g01 = (p00*ee + p01*ed) * (p10*ed + p11*ee);
				float g10 = (p00*ed + p01*ee) * (p10*ee + p11*ed);
				float sum = g01 + g10;
				g01 /= sum;
				g10 /= sum;

				if (abs(g01-g10) > threshold) {
					if (g01>g10) {
						GRvar_genotypes[vr][r].pha = 1;
						GRvar_genotypes[vr][r].al0 = 0;
						GRvar_genotypes[vr][r].al1 = 1;
					} else {
						GRvar_genotypes[vr][r].pha = 1;
						GRvar_genotypes[vr][r].al0 = 1;
						GRvar_genotypes[vr][r].al1 = 0;
					}
					n_phased ++;
				}
				n_total ++;
			}

			if (GRvar_genotypes[vr][r].mis) {
				float p00 = 1.0f - GRvar_genotypes[vr][r].ph0;
				float p01 = GRvar_genotypes[vr][r].ph0;
				float p10 = 1.0f - GRvar_genotypes[vr][r].ph1;
				float p11 = GRvar_genotypes[vr][r].ph1;
				gprobs[0] = (p00*ee + p01*ed) * (p10*ee + p11*ed);
				gprobs[1] = (p00*ee + p01*ed) * (p10*ed + p11*ee);
				gprobs[2] = (p00*ed + p01*ee) * (p10*ee + p11*ed);
				gprobs[3] = (p00*ed + p01*ee) * (p10*ed + p11*ee);
				int maxg = alg.imax(gprobs);
				GRvar_genotypes[vr][r].pha = 1;
				switch (maxg) {
				case 0:	GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 0; break;
				case 1:	GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 1; break;
				case 2:	GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 0; break;
				case 3:	GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 1; break;
				}
				n_phased ++;
				n_total ++;
			}
		}
	}
}

void genotype_set::allocate(variant_map & V, unsigned int _n_samples, unsigned int _n_scaffold_variants, unsigned int _n_rare_variants) {
	tac.clock();

	n_scaffold_variants = _n_scaffold_variants;
	n_rare_variants = _n_rare_variants;
	n_samples = _n_samples;

	if (n_rare_variants > 0) {
		GRvar_genotypes = vector < vector < rare_genotype > > (n_rare_variants);
		GRind_genotypes = vector < vector < rare_genotype > > (n_samples);
		MAP_R2S = vector < unsigned int > (n_rare_variants);
		major_alleles = vector < bool > (n_rare_variants, false);
		for (int r = 0 ; r < V.sizeRare() ; r ++) major_alleles[r] = !V.vec_rare[r]->minor;
	}

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

void genotype_set::transpose() {
	tac.clock();

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

void genotype_set::mapUnphasedOntoScaffold(int ind, vector < vector < unsigned int > > & map) {
	map.clear();
	map = vector < vector < unsigned int > > (n_scaffold_variants+1, vector < unsigned int > ());

	//Rare
	for (int vr = 0 ; vr < GRind_genotypes[ind].size() ; vr ++) {
		unsigned int idx = GRind_genotypes[ind][vr].idx;
		bool phased = GRind_genotypes[ind][vr].pha;
		bool missing = GRind_genotypes[ind][vr].mis;
		bool het = GRind_genotypes[ind][vr].het;
		if (!phased) {
			if (missing) map[MAP_R2S[idx]].push_back(idx);
			else if (het) map[MAP_R2S[idx]].push_back(idx);
		}
	}
}

void genotype_set::scaffoldTrio(int ikid, int ifather, int imother, vector < unsigned int > &counts) {
	for (int vr_kid = 0, vr_fat = 0, vr_mot = 0 ; vr_kid < GRind_genotypes[ikid].size() ; vr_kid ++) {

		unsigned int vr =  GRind_genotypes[ikid][vr_kid].idx;

		while (GRind_genotypes[ifather][vr_fat].idx < vr) vr_fat ++;
		while (GRind_genotypes[imother][vr_mot].idx < vr) vr_mot ++;

		if (GRind_genotypes[ikid][vr_kid].het) {
			bool father_is_hom = (GRind_genotypes[ifather][vr_fat].idx != vr) || ((GRind_genotypes[ifather][vr_fat].het + GRind_genotypes[ifather][vr_fat].mis)==0);
			bool mother_is_hom = (GRind_genotypes[imother][vr_mot].idx != vr) || ((GRind_genotypes[imother][vr_mot].het + GRind_genotypes[imother][vr_mot].mis)==0);

			if (father_is_hom && mother_is_hom) {
				bool fath0 = (GRind_genotypes[ifather][vr_fat].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				bool moth0 = (GRind_genotypes[imother][vr_mot].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				if (fath0 != moth0) {
					GRind_genotypes[ikid][vr_kid].pha = 1;
					GRind_genotypes[ikid][vr_kid].al0 = fath0;
					GRind_genotypes[ikid][vr_kid].al1 = moth0;
				} else counts[0]++;
			} else if (father_is_hom) {
				bool fath0 = (GRind_genotypes[ifather][vr_fat].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].al0 = fath0;
				GRind_genotypes[ikid][vr_kid].al1 = 1-fath0;
			} else if (mother_is_hom) {
				bool moth0 = (GRind_genotypes[imother][vr_mot].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].al0 = 1-moth0;
				GRind_genotypes[ikid][vr_kid].al1 = moth0;
			} else counts[3]++;
			counts[1] ++;
		} else if (GRind_genotypes[ikid][vr_kid].mis) {
			bool father_is_hom = (GRind_genotypes[ifather][vr_fat].idx != vr) || ((GRind_genotypes[ifather][vr_fat].het + GRind_genotypes[ifather][vr_fat].mis)==0);
			bool mother_is_hom = (GRind_genotypes[imother][vr_mot].idx != vr) || ((GRind_genotypes[imother][vr_mot].het + GRind_genotypes[imother][vr_mot].mis)==0);
			if (father_is_hom && mother_is_hom) {
				bool fath0 = (GRind_genotypes[ifather][vr_fat].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				bool moth0 = (GRind_genotypes[imother][vr_mot].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].al0 = fath0;
				GRind_genotypes[ikid][vr_kid].al1 = moth0;
			}
		}
	}
}

void genotype_set::scaffoldDuoMother(int ikid, int imother, vector < unsigned int > &counts) {
	for (int vr_kid = 0, vr_fat = 0, vr_mot = 0 ; vr_kid < GRind_genotypes[ikid].size() ; vr_kid ++) {
		unsigned int vr =  GRind_genotypes[ikid][vr_kid].idx;
		while (GRind_genotypes[imother][vr_mot].idx < vr) vr_mot ++;
		if (GRind_genotypes[ikid][vr_kid].het) {
			bool mother_is_hom = (GRind_genotypes[imother][vr_mot].idx != vr) || ((GRind_genotypes[imother][vr_mot].het + GRind_genotypes[imother][vr_mot].mis)==0);
			if (mother_is_hom) {
				bool moth0 = (GRind_genotypes[imother][vr_mot].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].al0 = 1-moth0;
				GRind_genotypes[ikid][vr_kid].al1 = moth0;
			} else counts[3]++;
			counts[1] ++;
		}
	}
}

void genotype_set::scaffoldDuoFather(int ikid, int ifather, vector < unsigned int > &counts) {
	for (int vr_kid = 0, vr_fat = 0 ; vr_kid < GRind_genotypes[ikid].size() ; vr_kid ++) {
		unsigned int vr =  GRind_genotypes[ikid][vr_kid].idx;
		while (GRind_genotypes[ifather][vr_fat].idx < vr) vr_fat ++;
		if (GRind_genotypes[ikid][vr_kid].het) {
			bool father_is_hom = (GRind_genotypes[ifather][vr_fat].idx != vr) || ((GRind_genotypes[ifather][vr_fat].het + GRind_genotypes[ifather][vr_fat].mis)==0);
			if (father_is_hom) {
				bool fath0 = (GRind_genotypes[ifather][vr_fat].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].al0 = fath0;
				GRind_genotypes[ikid][vr_kid].al1 = 1-fath0;
			} else counts[3]++;
			counts[1] ++;
		}
	}
}



//counts[0] : # observed mendel errors
//counts[1] : # possible mendel errors
//counts[2] : # hets being scaffolded
//counts[2] : # hets not being scaffolded
void genotype_set::scaffoldUsingPedigrees(pedigree_reader & pr) {
	tac.clock();
	vector < unsigned int > counts = vector < unsigned int >(4, 0);

	// Build map
	map < string, unsigned int > mapG;
	for (int i = 0 ; i < n_samples ; i ++) mapG.insert(pair < string, unsigned int > (names[i], i));

	//Mapping samples
	unsigned int ntrios = 0, nduos = 0, nmendels = 0;
	map < string, unsigned int > :: iterator itK, itM, itF;
	for (int i = 0 ; i < pr.kids.size() ; i ++) {
		itK = mapG.find(pr.kids[i]);
		itF = mapG.find(pr.fathers[i]);
		itM = mapG.find(pr.mothers[i]);
		int gkid = (itK != mapG.end())?itK->second : -1;
		int gfather = (itF != mapG.end())?itF->second : -1;
		int gmother = (itM != mapG.end())?itM->second : -1;
		if (gkid>=0) {
			if (gfather>=0 && gmother>=0) {
				scaffoldTrio(gkid, gfather, gmother, counts);
				ntrios++;
			} else if (gfather>=0) {
				scaffoldDuoFather(gkid, gfather, counts);
				nduos++;
			} else if (gmother>=0) {
				scaffoldDuoMother(gkid, gmother, counts);
				nduos++;
			}
		}
	}

	//Transpose


	//Verbose
	vrb.bullet("PED mapping (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("#trios = " + stb.str(ntrios) + " / #duos = " + stb.str(nduos));
	vrb.bullet2("%mendel_errors at kid_hets = " + stb.str(counts[0] *100.0 / counts[1], 2) + "% (n=" + stb.str(counts[0]) + ")");
	vrb.bullet2("%hets_phased = " + stb.str(counts[2]*100.0 / (counts[2]+counts[3]), 2) + "% (n=" + stb.str(counts[2]) + ")");
}


