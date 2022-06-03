/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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
#ifndef _GENOTYPE_SET_H
#define _GENOTYPE_SET_H

#include <utils/otools.h>

#include <containers/bitmatrix.h>
#include <containers/variant_map.h>

#include <io/pedigree_reader.h>


class rare_genotype {
public:

	static float ee;
	static float ed;

	unsigned int idx : 27;
	unsigned int het : 1;
	unsigned int mis : 1;

	unsigned int pha : 1;
	unsigned int al0 : 1;
	unsigned int al1 : 1;

	float ph0, ph1;

	rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
		ph0 = ph1 = 0.0f;
	}

	rare_genotype(unsigned int _idx, bool _het, bool _mis, bool _al0, bool _al1) {
		idx = _idx; het = _het; mis = _mis; al0 = _al0; al1 = _al1;
		pha = (!mis && (al0==al1));
		ph0 = ph1 = 0.0f;
	}

	~rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
		ph0 = ph1 = 0.0f;
	}

	bool operator < (const rare_genotype & rg) const {
		return idx < rg.idx;
	}

	void phase(unsigned int g) {
		pha = 1;
		switch (g) {
		case 0:	al0 = 0; al1 = 0; break;
		case 1:	al0 = 0; al1 = 1; break;
		case 2:	al0 = 1; al1 = 0; break;
		case 3:	al0 = 1; al1 = 1; break;
		}
	}

	void randomize() {
		if (het) {
			bool flip = rng.flipCoin();
			al0 = flip?true:false;
			al1 = flip?false:true;
		}

		if (mis) {
			switch (rng.getInt(4)) {
			case 0:	al0 = false; al1 = false; break;
			case 1:	al0 = false; al1 = true; break;
			case 2:	al0 = true; al1 = false; break;
			case 3:	al0 = true; al1 = true; break;
			}
		}
	}

	void computeProbs(vector < float > & gprobs) {
		assert(gprobs.size() == 4);

		float p00 = 1.0f - ph0;
		float p01 = ph0;
		float p10 = 1.0f - ph1;
		float p11 = ph1;

		if (het) {
			gprobs[0] = 0.0f;
			gprobs[1] = (p00*ee + p01*ed) * (p10*ed + p11*ee);
			gprobs[2] = (p00*ed + p01*ee) * (p10*ee + p11*ed);
			gprobs[3] = 0.0f;
		} else if (mis) {
			gprobs[0] = (p00*ee + p01*ed) * (p10*ee + p11*ed);
			gprobs[1] = (p00*ee + p01*ed) * (p10*ed + p11*ee);
			gprobs[2] = (p00*ed + p01*ee) * (p10*ee + p11*ed);
			gprobs[3] = (p00*ed + p01*ee) * (p10*ed + p11*ee);
		}

		float scale = 1.0f / (gprobs[0] + gprobs[1] + gprobs[2] + gprobs[3]);

		gprobs[0] *= scale;
		gprobs[1] *= scale;
		gprobs[2] *= scale;
		gprobs[3] *= scale;
	}
};

class genotype_set {
public:

	//Counts
	unsigned int n_scaffold_variants;			//#variants rare to be phased
	unsigned int n_rare_variants;				//#variants rare to be phased
	unsigned int n_samples;						//#samples

	//Sample IDs
	vector < string > names;

	//Mapping on scaffold
	vector < unsigned int > MAP_R2S;

	//Genotypes at rare unphased variants
	vector < bool > major_alleles;
	vector < vector < rare_genotype > > GRvar_genotypes;
	vector < vector < rare_genotype > > GRind_genotypes;

	//
	genotype_set();
	~genotype_set();
	void clear();
	void allocate(variant_map &, unsigned int, unsigned int , unsigned int);
	void transpose();
	void mapUnphasedOntoScaffold(int ind, vector < vector < unsigned int > > & map);
	void impute(unsigned int, unsigned int, vector < float > &, vector < float > &, vector < unsigned int > &);

	//
	void pushRareMissing(unsigned int vr, unsigned int i, bool major);
	void pushRareHet(unsigned int vr, unsigned int i);
	void pushRareHom(unsigned int vr, unsigned int i, bool major);
	void imputeMonomorphic();
	void randomizeSingleton();
	void phase(unsigned long int & n_phased, unsigned long int & n_total, float);

	void scaffoldTrio(int ikid, int ifather, int imother, vector < unsigned int > &counts);
	void scaffoldDuoMother(int ikid, int imother, vector < unsigned int > &counts);
	void scaffoldDuoFather(int ikid, int ifather, vector < unsigned int > &counts);
	void scaffoldUsingPedigrees(pedigree_reader & pr);
};

inline
void genotype_set::pushRareMissing(unsigned int vr, unsigned int i, bool major) {
	GRvar_genotypes[vr].emplace_back(i, 0, 1, major, major);
}

inline
void genotype_set::pushRareHet(unsigned int vr, unsigned int i) {
	GRvar_genotypes[vr].emplace_back(i, 1, 0, 0, 1);
}

inline
void genotype_set::pushRareHom(unsigned int vr, unsigned int i, bool major) {
	GRvar_genotypes[vr].emplace_back(i, 0, 0, !major, !major);
}

#endif
