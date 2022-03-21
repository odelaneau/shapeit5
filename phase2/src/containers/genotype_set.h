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


class rare_genotype {
public:
	unsigned int idx : 28;
	unsigned int het : 1;
	unsigned int mis : 1;
	unsigned int al0 : 1;
	unsigned int al1 : 1;

	rare_genotype() { idx = het = mis = al0 = al1 = 0; }

	rare_genotype(unsigned int _idx, bool _het, bool _mis, bool _al0, bool _al1) {
		idx = _idx; het = _het; mis = _mis; al0 = _al0; al1 = _al1;
	}

	~rare_genotype() {
		idx = het = mis = al0 = al1 = 0;
	}

	bool operator < (const rare_genotype & rg) const {
		return idx < rg.idx;
	}
};

class genotype_set {
public:

	//Counts
	unsigned int n_scaffold_variants;				//#variants rare to be phased
	unsigned int n_rare_variants;				//#variants rare to be phased
	unsigned int n_common_variants;				//#variants common to be phased (e.g. indels)
	unsigned int n_samples;						//#samples

	//Sample IDs
	vector < string > names;

	//Mapping on scaffold
	vector < unsigned int > MAPC_vs_left;
	vector < unsigned int > MAPC_vs_right;
	vector < unsigned int > MAPR_vs_left;
	vector < unsigned int > MAPR_vs_right;

	//Genotypes at common unphased variants
	bitmatrix GCvar_alleles;
	bitmatrix GCind_alleles;
	bitmatrix GCvar_missing;
	bitmatrix GCind_missing;

	//Genotypes at rare unphased variants
	vector < vector < rare_genotype > > GRvar_genotypes;
	vector < vector < rare_genotype > > GRind_genotypes;

	//
	genotype_set();
	~genotype_set();
	void clear();
	void allocate(variant_map &, unsigned int, unsigned int , unsigned int, unsigned int);
	void transpose();
	void mapUnphasedOntoScaffold(int ind, vector < bool > & map);

	//
	void setCommonMissing(unsigned int vc, unsigned int i);
	void setCommonGenotype(unsigned int vc, unsigned int i, unsigned int g);
	void pushRareMissing(unsigned int vr, unsigned int i, bool major);
	void pushRareHet(unsigned int vr, unsigned int i);
	void pushRareHom(unsigned int vr, unsigned int i, bool major);
};

inline
void genotype_set::setCommonMissing(unsigned int vc, unsigned int i) {
	GCvar_missing.set(vc, i, true);
	GCvar_alleles.set(vc, 2*i+0, false);
	GCvar_alleles.set(vc, 2*i+1, false);
}

inline
void genotype_set::setCommonGenotype(unsigned int vc, unsigned int i, unsigned int g) {
	GCvar_missing.set(vc, i, false);
	if (g == 0) {
		GCvar_alleles.set(vc, 2*i+0, false);
		GCvar_alleles.set(vc, 2*i+1, false);
	} else if (g == 2) {
		GCvar_alleles.set(vc, 2*i+0, true);
		GCvar_alleles.set(vc, 2*i+1, true);
	} else {
		GCvar_alleles.set(vc, 2*i+0, false);
		GCvar_alleles.set(vc, 2*i+1, true);
	}
}

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
