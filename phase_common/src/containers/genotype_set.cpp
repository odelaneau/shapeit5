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

#include <containers/genotype_set.h>

using namespace std;

genotype_set::genotype_set() {
	n_site = 0;
	n_ind = 0;
}

genotype_set::~genotype_set() {
	for (int i = 0 ; i< vecG.size() ; i ++) delete vecG[i];
	vecG.clear();
	n_site = 0;
	n_ind = 0;
}

void genotype_set::allocate(unsigned long n_main_samples, unsigned long n_variants) {
	vecG = vector < genotype * > (n_main_samples);
	for (unsigned int i = 0 ; i < n_main_samples ; i ++) {
		vecG[i] = new genotype (i);
		vecG[i]->n_variants = n_variants;
		vecG[i]->Variants = vector < unsigned char > (DIV2(n_variants) + MOD2(n_variants), 0);
	}
	n_ind = n_main_samples;
	n_site = n_variants;
}

void genotype_set::imputeMonomorphic(variant_map & V) {
	tac.clock();
	unsigned int n_imputed_genotypes = 0;
	for (unsigned int v = 0 ; v < V.size() ; v ++) {
		if (V.vec_pos[v]->isMonomorphic()) {
			bool uallele = (V.vec_pos[v]->cref)?false:true;
			for (unsigned int i = 0 ; i < vecG.size() ; i ++) {
				VAR_SET_HOM(MOD2(v), vecG[i]->Variants[DIV2(v)]);
				uallele?VAR_SET_HAP0(MOD2(v), vecG[i]->Variants[DIV2(v)]):VAR_CLR_HAP0(MOD2(v), vecG[i]->Variants[DIV2(v)]);
				uallele?VAR_SET_HAP1(MOD2(v), vecG[i]->Variants[DIV2(v)]):VAR_CLR_HAP1(MOD2(v), vecG[i]->Variants[DIV2(v)]);
				n_imputed_genotypes++;
			}
			if (uallele) V.vec_pos[v]->cref = 0;
			else V.vec_pos[v]->calt = 0;
			V.vec_pos[v]->cmis = 0;
		}
	}
	vrb.bullet("Impute monomorphic [n=" + stb.str(n_imputed_genotypes) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

unsigned int genotype_set::largestNumberOfTransitions() {
	unsigned int maxT = 0;
	for (int i = 0 ; i < n_ind ; i ++) {
		unsigned int nTrans = vecG[i]->n_transitions;
		if (nTrans > maxT) maxT = nTrans;
	}
	return maxT;
}

unsigned int genotype_set::largestNumberOfMissings() {
	unsigned int maxM = 0;
	for (int i = 0 ; i < n_ind ; i ++) {
		unsigned int nMis = vecG[i]->n_missing * HAP_NUMBER;
		if (nMis> maxM) maxM = nMis;
	}
	return maxM;
}

unsigned long genotype_set::numberOfSegments() {
	unsigned long size = 0;
	for (int i = 0 ; i < n_ind ; i ++) size += vecG[i]->n_segments;
	return size;
}

void genotype_set::solve() {
	tac.clock();
	for (int i = 0 ; i < vecG.size() ; i ++) vecG[i]->solve();
	vrb.bullet("HAP solving (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

//counts[0] : # observed mendel errors
//counts[1] : # possible mendel errors
//counts[2] : # hets being scaffolded
//counts[3] : # hets not being scaffolded
void genotype_set::scaffoldUsingPedigrees(pedigree_reader & pr) {
	tac.clock();
	vector < unsigned int > counts = vector < unsigned int >(4, 0);

	// Build map
	map < string, genotype * > mapG;
	for (int i = 0 ; i < n_ind ; i ++) mapG.insert(pair < string, genotype * > (vecG[i]->name, vecG[i]));

	//Mapping samples
	unsigned int ntrios = 0, nduos = 0, nmendels = 0;
	map < string, genotype * > :: iterator itK, itM, itF;
	for (int i = 0 ; i < pr.kids.size() ; i ++) {
		itK = mapG.find(pr.kids[i]);
		itF = mapG.find(pr.fathers[i]);
		itM = mapG.find(pr.mothers[i]);
		genotype * gkid = (itK != mapG.end())?itK->second : NULL;
		genotype * gfather = (itF != mapG.end())?itF->second : NULL;
		genotype * gmother = (itM != mapG.end())?itM->second : NULL;
		if (gkid) {
			if (gfather && gmother) {
				gkid->scaffoldTrio(gfather, gmother, counts);
				ntrios++;
			} else if (gfather) {
				gkid->scaffoldDuoFather(gfather, counts);
				nduos++;
			} else if (gmother) {
				gkid->scaffoldDuoMother(gmother, counts);
				nduos++;
			}
		}
	}

	//Verbose
	vrb.bullet("PED mapping (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("#trios = " + stb.str(ntrios) + " / #duos = " + stb.str(nduos));
	vrb.bullet2("%mendel_errors = " + stb.str(counts[0] *100.0 / counts[1], 2) + "% (n=" + stb.str(counts[0]) + ")");
	vrb.bullet2("%hets_phased = " + stb.str(counts[2]*100.0 / (counts[2]+counts[3]), 2) + "% (n=" + stb.str(counts[2]) + ")");
}

