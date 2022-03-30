/*******************************************************************************
 * Copyright (C) 2018-2022 Olivier Delaneau, University of Lausanne
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
#ifndef _GENOTYPE_H
#define _GENOTYPE_H

#include <utils/otools.h>

//Macros for packing/unpacking diplotypes
#define DIP_GET(dip,idx)	(((dip)>>(idx))&1UL)
#define DIP_SET(dip,idx)	((dip)|=(1UL<<(idx)))
#define DIP_HAP0(idx)		((idx)>>3)
#define DIP_HAP1(idx)		((idx)&7)

//Macros for packing/unpacking haplotypes
#define HAP_GET(hap, idx)	(((hap)>>(idx))&1U)
#define HAP_SET(hap, idx)	((hap)|=(1U<<(idx)))

//Macros for packing/unpacking variants
#define VAR_GET_HOM(e,v)	((((v)>>((e)<<2))&3)==0)
#define VAR_GET_MIS(e,v)	((((v)>>((e)<<2))&3)==1)
#define VAR_GET_HET(e,v)	((((v)>>((e)<<2))&3)==2)
#define VAR_GET_SCA(e,v)	((((v)>>((e)<<2))&3)==3)
//#define VAR_GET_AMB(e,v)	((((v)>>((e)<<2))&3)!=0)
#define VAR_GET_AMB(e,v)	((((v)>>((e)<<2))&3)>1)
#define VAR_SET_HOM(e,v)	((e)?((v)&=207):((v)&=252))
#define VAR_SET_MIS(e,v)	((v)|=(1<<((e)<<2)))
#define VAR_SET_HET(e,v)	((v)|=(2<<((e)<<2)))
#define VAR_SET_SCA(e,v)	((v)|=(3<<((e)<<2)))

#define VAR_GET_HAP0(e,v)	(((v)&(4<<((e)<<2)))!=0)
#define VAR_SET_HAP0(e,v)	((v)|=(4<<((e)<<2)))
#define VAR_CLR_HAP0(e,v)	((e)?((v)&=191):((v)&=251))
#define VAR_GET_HAP1(e,v)	(((v)&(8<<((e)<<2)))!=0)
#define VAR_SET_HAP1(e,v)	((v)|=(8<<((e)<<2)))
#define VAR_CLR_HAP1(e,v)	((e)?((v)&=127):((v)&=247))

class genotype {
public:
	// INTERNAL DATA
	string name;
	unsigned int index;						// Index in containers
	unsigned int n_segments;				// Number of segments
	unsigned int n_variants;				// Number of variants	(to iterate over Variants)
	unsigned int n_ambiguous;				// Number of ambiguous variants
	unsigned int n_missing;					// Number of missing
	unsigned int n_transitions;				// Number of transitions
	unsigned int n_stored_transitionProbs;	// Number of transition probabilities stored in memory
	unsigned int n_storage_events;			// Number of storage having been done
	unsigned char curr_dipcodes [64];		// List of diplotypes in a given segment (buffer style variable)

	// VARIANT / HAPLOTYPE / DIPLOTYPE DATA
	vector < unsigned char > Variants;		// 0.5 byte per variant
	vector < unsigned char > Ambiguous;		// 1 byte per ambiguous variant
	vector < unsigned long > Diplotypes;	// 8 bytes per segment
	vector < unsigned short > Lengths;		// 2 bytes per segment

	//PHASE PROBS
	vector < bool > ProbMask;
	vector < float > ProbStored;
	vector < float > ProbMissing;

	//METHODS
	genotype(unsigned int);
	~genotype();
	void free();
	void make(vector < unsigned char > &, vector < float > &);
	void make(vector < unsigned char > &);
	void build();
	void sample(vector < double > &, vector < float > &);
	void sampleForward(vector < double > &, vector < float > &);
	void sampleBackward(vector < double > &, vector < float > &);
	void solve();
	void mapMerges(vector < double > &, double , vector < bool > &);
	void performMerges(vector < double > &, vector < bool > &);
	void store(vector < double > &, vector < float > &);
	void scaffoldTrio(genotype *, genotype *, vector < unsigned int > &);
	void scaffoldDuoFather(genotype *, vector < unsigned int > &);
	void scaffoldDuoMother(genotype *, vector < unsigned int > &);

	//INLINES
	unsigned int countDiplotypes(unsigned long);
	void makeDiplotypes(unsigned long);
	unsigned int countTransitions();
};

inline
unsigned int genotype::countDiplotypes(unsigned long _dip) {
	unsigned int c = 0;
	for (unsigned long dip = _dip; dip; c++) dip &= dip - 1;
	return c;
}

inline
void genotype::makeDiplotypes(unsigned long _dip) {
	for (unsigned int d = 0, i = 0 ; d < 64 ; ++d) if (DIP_GET(_dip, d)) curr_dipcodes[i++] = d;
}

inline
unsigned int genotype::countTransitions() {
	unsigned int prev_dipcount = 1, c = 0;
	for (unsigned int s = 0 ; s < n_segments ; s++) {
		unsigned int curr_dipcount = countDiplotypes(Diplotypes[s]);
		c+= prev_dipcount * curr_dipcount;
		prev_dipcount = curr_dipcount;
	}
	return c;
}

#endif
