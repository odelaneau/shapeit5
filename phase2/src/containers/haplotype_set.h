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
#ifndef _HAPLOTYPE_SET_H
#define _HAPLOTYPE_SET_H

#include <utils/otools.h>

#include <containers/bitmatrix.h>
#include <containers/genotype_set.h>
#include <containers/variant_map.h>

class haplotype_set {
public:
	//
	unsigned int n_common_variants;
	unsigned int n_rare_variants;
	unsigned int n_total_variants;
	unsigned int n_total_haplotypes;
	unsigned int n_target_haplotypes;
	unsigned int n_reference_haplotypes;
	unsigned int n_total_samples;
	unsigned int n_target_samples;
	unsigned int n_reference_samples;

	//Haplotype Data [COMMONS]
	bitmatrix Hhap;					// Bit matrix of haplotypes (haplotype first). Transposed version of H_opt_var.
	bitmatrix Hvar;					// Bit matrix of haplotypes (variant first). Transposed version of H_opt_hap.

	//Haplotype Data [RARES]
	vector < vector < unsigned int > > Shap;
	vector < vector < unsigned int > > Svar;
	vector < bool > flag_common;

	//PBWT parameters
	double pbwt_modulo;				// Modulo used to store PBWT indexes (--pbwt-modulo)
	unsigned long pbwt_depth;		// #neighbours in the PBWT to use for conditioning (--pbwt-depth)
	unsigned long pbwt_mac;			// Minor Allele Count to consider in PBWT pass (--pbwt-mac)
	double pbwt_mdr;				// Missinga Data Rate to consider in PBWT pass (--pbwt-mdr)
	unsigned int nthreads;			// Number of threads (--thread)

	//PBWT indexes & arrays
	unsigned long pbwt_nstored;		//#variants with PBWT indexes stored
	vector < double > pbwt_cm;		//Variants at which PBWT is evaluated
	vector < int > pbwt_grp;		//Variant groups based on cm positions
	vector < int > pbwt_evaluated;	//Variants at which PBWT is evaluated
	vector < int > pbwt_stored;		//Variants at which PBWT is stored
	vector < int > pbwt_parray;		//PBWT prefix array
	vector < int > pbwt_darray;		//PBWT divergence array
	//vector < int > pbwt_neighbours; //Closest neighbours
	vector < vector < unsigned int > > pbwt_neighbours;

	//PBWT IBD2 protect
	vector < vector < unsigned int > > bannedPairs;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();
	void clear();

	//
	void allocate(unsigned int,unsigned int , unsigned int , unsigned int, vector < bool > &);

	//PBWT routines
	void parametrizePBWT(int, double, int, double, int);
	void initializePBWTmapping(variant_map &);
	void updatePBWTmapping();
	void allocatePBWTarrays();
	void selectPBWTarrays();
	//void transposePBWTarrays();

	//IBD2 routines
	void mergeIBD2constraints();
	bool checkIBD2matching(int, int);

	//Haplotype routines
	void updateHaplotypes(genotype_set & G, bool first_time = false);
	void transposeHaplotypes_H2V(bool full);
	void transposeHaplotypes_V2H(bool full);
};


inline
bool haplotype_set::checkIBD2matching(int mh, int ch) {
	int mi = min(mh/2,ch/2);
	int ci = max(mh/2,ch/2);
	// Prevents self copying, who knows ?
	if (mi == ci) return false;
	// Prevents copying for IBD2 individuals
	for (int i = 0 ; i < bannedPairs[mi].size() && bannedPairs[mi][i] <= ci; i ++) if (bannedPairs[mi][i] == ci) return false;
	return true;
}

#endif
