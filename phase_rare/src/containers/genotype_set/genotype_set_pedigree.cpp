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

void genotype_set::phaseTrio(int ikid, int ifather, int imother, vector < unsigned int > &counts) {
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
					GRind_genotypes[ikid][vr_kid].prob = 1.0f;
					GRind_genotypes[ikid][vr_kid].al0 = fath0;
					GRind_genotypes[ikid][vr_kid].al1 = moth0;
				} else counts[0]++;
			} else if (father_is_hom) {
				bool fath0 = (GRind_genotypes[ifather][vr_fat].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].prob = 1.0f;
				GRind_genotypes[ikid][vr_kid].al0 = fath0;
				GRind_genotypes[ikid][vr_kid].al1 = 1-fath0;
			} else if (mother_is_hom) {
				bool moth0 = (GRind_genotypes[imother][vr_mot].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].prob = 1.0f;
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
				GRind_genotypes[ikid][vr_kid].prob = 1.0f;
				GRind_genotypes[ikid][vr_kid].al0 = fath0;
				GRind_genotypes[ikid][vr_kid].al1 = moth0;
			}
		}
	}
}

void genotype_set::phaseDuoMother(int ikid, int imother, vector < unsigned int > &counts) {
	for (int vr_kid = 0, vr_fat = 0, vr_mot = 0 ; vr_kid < GRind_genotypes[ikid].size() ; vr_kid ++) {
		unsigned int vr =  GRind_genotypes[ikid][vr_kid].idx;
		while (GRind_genotypes[imother][vr_mot].idx < vr) vr_mot ++;
		if (GRind_genotypes[ikid][vr_kid].het) {
			bool mother_is_hom = (GRind_genotypes[imother][vr_mot].idx != vr) || ((GRind_genotypes[imother][vr_mot].het + GRind_genotypes[imother][vr_mot].mis)==0);
			if (mother_is_hom) {
				bool moth0 = (GRind_genotypes[imother][vr_mot].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].prob = 1.0f;
				GRind_genotypes[ikid][vr_kid].al0 = 1-moth0;
				GRind_genotypes[ikid][vr_kid].al1 = moth0;
			} else counts[3]++;
			counts[1] ++;
		}
	}
}

void genotype_set::phaseDuoFather(int ikid, int ifather, vector < unsigned int > &counts) {
	for (int vr_kid = 0, vr_fat = 0 ; vr_kid < GRind_genotypes[ikid].size() ; vr_kid ++) {
		unsigned int vr =  GRind_genotypes[ikid][vr_kid].idx;
		while (GRind_genotypes[ifather][vr_fat].idx < vr) vr_fat ++;
		if (GRind_genotypes[ikid][vr_kid].het) {
			bool father_is_hom = (GRind_genotypes[ifather][vr_fat].idx != vr) || ((GRind_genotypes[ifather][vr_fat].het + GRind_genotypes[ifather][vr_fat].mis)==0);
			if (father_is_hom) {
				bool fath0 = (GRind_genotypes[ifather][vr_fat].idx != vr)?major_alleles[vr]:!major_alleles[vr];
				GRind_genotypes[ikid][vr_kid].pha = 1;
				GRind_genotypes[ikid][vr_kid].prob = 1.0f;
				GRind_genotypes[ikid][vr_kid].al0 = fath0;
				GRind_genotypes[ikid][vr_kid].al1 = 1-fath0;
			} else counts[3]++;
			counts[1] ++;
		}
	}
}

//counts[0] : # observed mendel errors
//counts[1] : # possible mendel errors
//counts[2] : # hets being phased
//counts[3] : # hets not being phased
void genotype_set::phaseUsingPedigrees(pedigree_reader & pr) {
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
				phaseTrio(gkid, gfather, gmother, counts);
				ntrios++;
			} else if (gfather>=0) {
				phaseDuoFather(gkid, gfather, counts);
				nduos++;
			} else if (gmother>=0) {
				phaseDuoMother(gkid, gmother, counts);
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


