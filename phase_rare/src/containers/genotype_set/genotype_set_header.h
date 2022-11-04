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

#ifndef _GENOTYPE_SET_H
#define _GENOTYPE_SET_H

#include <utils/otools.h>

#include <containers/bitmatrix.h>
#include <objects/hmm_parameters.h>
#include <objects/rare_genotype.h>
#include <io/pedigree_reader.h>

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
	std::vector < std::string > names;

	//Mapping on scaffold
	std::vector < unsigned int > MAP_R2S;

	//Genotypes at rare unphased variants
	std::vector < bool > major_alleles;
	std::vector < std::vector < rare_genotype > > GRvar_genotypes;
	std::vector < std::vector < rare_genotype > > GRind_genotypes;

	//
	genotype_set();
	~genotype_set();

	//BASIC
	void clear();
	void allocate(variant_map &, unsigned int, unsigned int , unsigned int);
	void mapUnphasedOntoScaffold(int ind, std::vector < std::vector < unsigned int > > & map);

	//TRANSPOSE
	void fillup_by_transpose_V2I();
	void merge_by_transpose_I2V();


	//PUSH
	void pushRareMissing(unsigned int vr, unsigned int i, bool major);
	void pushRareHet(unsigned int vr, unsigned int i);
	void pushRareHom(unsigned int vr, unsigned int i, bool major);
	int pushRare(unsigned int vr, unsigned int v);

	//IMPUTE
	void imputeMonomorphic();
	void phaseLiAndStephens(unsigned int, unsigned int, aligned_vector32 < float > &, aligned_vector32 < float > &, std::vector < unsigned int > &, float);
	void phaseCoalescentViterbi(unsigned int, std::vector < int > &, std::vector < int > &, hmm_parameters &);

	//TRIOS [UNTESTED]
	void phaseTrio(int ikid, int ifather, int imother, std::vector < unsigned int > &counts);
	void phaseDuoMother(int ikid, int imother, std::vector < unsigned int > &counts);
	void phaseDuoFather(int ikid, int ifather, std::vector < unsigned int > &counts);
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

inline
int genotype_set::pushRare(unsigned int vr, unsigned int value) {
	GRvar_genotypes[vr].emplace_back(value);
	if (GRvar_genotypes[vr].back().mis) return 3;
	else if (GRvar_genotypes[vr].back().het) return 1;
	else return 2*GRvar_genotypes[vr].back().al0;
}

#endif
