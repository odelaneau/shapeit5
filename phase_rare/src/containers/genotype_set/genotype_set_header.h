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
#include <objects/hmm_parameters.h>

#include <io/pedigree_reader.h>


class rare_genotype {
public:

	static float ee;
	static float ed;

	unsigned int idx : 27;	//Index of the non Major/Major genotype
	unsigned int het : 1;	//Is it het?
	unsigned int mis : 1;	//Is it missing?	If none of the two, it is then Minor/Minor
	unsigned int al0 : 1;	//What's the allele on haplotype 0?
	unsigned int al1 : 1;	//What's the allele on haplotype 1?
	unsigned int pha : 1;	//Is the genotype phased?
	float prob;				//Probability

	rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
		prob = -1.0f;
	}

	rare_genotype(unsigned int _idx, bool _het, bool _mis, bool _al0, bool _al1) {
		idx = _idx; het = _het; mis = _mis; al0 = _al0; al1 = _al1;
		pha = (!mis && (al0==al1));

		if (al0 != al1) {
			if (rng.flipCoin()) {
				al0 = 0;
				al1 = 1;
			} else {
				al0 = 1;
				al1 = 0;
			}
		}

		prob = pha?1.0f:-1.0f;
	}

	~rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
		prob = -1.0f;
	}

	bool operator < (const rare_genotype & rg) const {
		return idx < rg.idx;
	}

	void phase(unsigned int g) {
		if (!pha) {
			switch (g) {
			case 0:	al0 = 0; al1 = 0; break;
			case 1:	al0 = 0; al1 = 1; break;
			case 2:	al0 = 1; al1 = 0; break;
			case 3:	al0 = 1; al1 = 1; break;
			}
		}
	}

	void randomize() {
		if (!pha) {
			if (het) phase(rng.getInt(2) + 1);
			if (mis) phase(rng.getInt(4));
		}
	}

	void phase(float prb0, float prb1) {
		if (!pha && (het || mis)) {
			vector < double > gprobs = vector < double > (4, 0.0f);


			float p01 = max(prb0, numeric_limits<float>::min());
			float p00 = max(1.0f - prb0, numeric_limits<float>::min());
			float p11 = max(prb1, numeric_limits<float>::min());
			float p10 = max(1.0f - prb1, numeric_limits<float>::min());

			if (het) {
				gprobs[0] = 0.0f;
				gprobs[1] = (p00*ee + p01*ed) * (p10*ed + p11*ee);
				gprobs[2] = (p00*ed + p01*ee) * (p10*ee + p11*ed);
				gprobs[3] = 0.0f;
			} else {
				gprobs[0] = (p00*ee + p01*ed) * (p10*ee + p11*ed);
				gprobs[1] = (p00*ee + p01*ed) * (p10*ed + p11*ee);
				gprobs[2] = (p00*ed + p01*ee) * (p10*ee + p11*ed);
				gprobs[3] = (p00*ed + p01*ee) * (p10*ed + p11*ee);
			}

			/*
			if (het) {
				gprobs[0] = 0.0f;
				gprobs[1] = p00 * p11;
				gprobs[2] = p01 * p10;
				gprobs[3] = 0.0f;
			} else {
				gprobs[0] = p00 * p10;
				gprobs[1] = p00 * p11;
				gprobs[2] = p01 * p10;
				gprobs[3] = p01 * p11;
			}
			*/

			//cout << stb.str(gprobs, 3) << endl;
			int maxg = alg.imax(gprobs);
			switch (maxg) {
			case 0:	al0 = 0; al1 = 0; break;
			case 1:	al0 = 0; al1 = 1; break;
			case 2:	al0 = 1; al1 = 0; break;
			case 3:	al0 = 1; al1 = 1; break;
			}
			prob = gprobs[maxg] / (gprobs[0] + gprobs[1] + gprobs[2] + gprobs[3]);
			//cout << stb.str(prob, 6) << endl;
		}
	}
};

class genotype_set {
public:

	//Counts
	unsigned int n_scaffold_variants;			//#variants rare to be phased
	unsigned int n_rare_variants;				//#variants rare to be phased
	unsigned int n_samples;						//#samples

	unsigned int nmiss_total;
	unsigned int nmiss_imputation;
	unsigned int nmiss_singleton;
	unsigned int nhets_total;
	unsigned int nhets_families;
	unsigned int nhets_imputation;
	unsigned int nhets_coalescent;

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

	//BASIC
	void clear();
	void allocate(variant_map &, unsigned int, unsigned int , unsigned int);
	void mapUnphasedOntoScaffold(int ind, vector < vector < unsigned int > > & map);

	//TRANSPOSE
	void fillup_by_transpose_V2I();
	void merge_by_transpose_I2V();


	//PUSH
	void pushRareMissing(unsigned int vr, unsigned int i, bool major);
	void pushRareHet(unsigned int vr, unsigned int i);
	void pushRareHom(unsigned int vr, unsigned int i, bool major);

	//IMPUTE
	void imputeMonomorphic();
	void phaseLiAndStephens(unsigned int, unsigned int, aligned_vector32 < float > &, aligned_vector32 < float > &, vector < unsigned int > &, float);
	void phaseCoalescentViterbi(unsigned int, vector < int > &, vector < int > &, hmm_parameters &);
	void phaseCoalescentViterbi2(unsigned int, double, double);

	//TRIOS [UNTESTED]
	void phaseTrio(int ikid, int ifather, int imother, vector < unsigned int > &counts);
	void phaseDuoMother(int ikid, int imother, vector < unsigned int > &counts);
	void phaseDuoFather(int ikid, int ifather, vector < unsigned int > &counts);
	void phaseUsingPedigrees(pedigree_reader & pr);
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
