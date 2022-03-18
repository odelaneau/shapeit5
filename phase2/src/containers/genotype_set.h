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

class genotype_set {
public:

	//Counts
	unsigned int n_scaffold_variants;			//#variants in scaffold
	unsigned int n_rare_variants;				//#variants rare to be phased
	unsigned int n_common_variants;				//#variants common to be phased (e.g. indels)
	unsigned int n_total_variants;				//#variants in total
	unsigned int n_samples;						//#samples

	//Genotypes at common unphased variants
	vector < vector < bool > > GCalleles;	//Alleles
	vector < vector < bool > > GCmissing;	//Missing?

	//Genotypes at rare unphased variants
	vector < vector < unsigned int > > GRindexes;	//Indexes of [ Maj/Min + Min/Min + ./.]
	vector < vector < bool > > GRhets;				//Maj/Min OR Min/Min?
	vector < vector < bool > > GRmissing;			//./.?
	vector < vector < bool > > GRalleles;			//

	//
	genotype_set();
	~genotype_set();
	void allocate(unsigned int,unsigned int , unsigned int , unsigned int, variant_map &);

	//
	void setCommonMissing(unsigned int vc, unsigned int i);
	void setCommonGenotype(unsigned int vc, unsigned int i, unsigned int g);
	void pushRareMissing(unsigned int vr, unsigned int i);
	void pushRareHet(unsigned int vr, unsigned int i);
	void pushRareHom(unsigned int vr, unsigned int i);

	//
	void getUnphasedIndexes(vector < unsigned int > &, vector < unsigned int > &, vector < unsigned int > &);
};
#endif

inline
void genotype_set::setCommonMissing(unsigned int vc, unsigned int i) {
	GCmissing[vc][i] = true;
}

inline
void genotype_set::setCommonGenotype(unsigned int vc, unsigned int i, unsigned int g) {
	if (g == 2) {
		GCalleles[vc][2*i+0] = true;
		GCalleles[vc][2*i+1] = true;
	} else if (g == 1) GCalleles[vc][2*i+1] = true;
}

inline
void genotype_set::pushRareMissing(unsigned int vr, unsigned int i) {
	GRindexes[vr].push_back(i);
	GRhets[vr].push_back(false);
	GRmissing[vr].push_back(true);
	GRalleles[vr].push_back(false);
	GRalleles[vr].push_back(false);
}

inline
void genotype_set::pushRareHet(unsigned int vr, unsigned int i) {
	GRindexes[vr].push_back(i);
	GRhets[vr].push_back(true);
	GRmissing[vr].push_back(false);
	GRalleles[vr].push_back(false);
	GRalleles[vr].push_back(false);
}

inline
void genotype_set::pushRareHom(unsigned int vr, unsigned int i) {
	GRindexes[vr].push_back(i);
	GRhets[vr].push_back(false);
	GRmissing[vr].push_back(false);
	GRalleles[vr].push_back(false);
	GRalleles[vr].push_back(false);
}
