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



#define MASK_INIT	0xFFFFFFFFFFFFFFFFUL
#define MASK_SCAF	0x00AA00AA00AA00AAUL
#define MASK_UNF0	0x55AA55AA55AA55AAUL
#define MASK_UNF1	0x3333CCCC3333CCCCUL
#define MASK_UNF2	0x0F0F0F0FF0F0F0F0UL


class genotype {
public:
	// INTERNAL DATA
	std::string name;
	unsigned int index;						// Index in containers
	unsigned int n_segments;				// Number of segments
	unsigned int n_variants;				// Number of variants	(to iterate over Variants)
	unsigned int n_ambiguous;				// Number of ambiguous variants
	unsigned int n_missing;					// Number of missing
	unsigned int n_transitions;				// Number of transitions
	unsigned int n_stored_transitionProbs;	// Number of transition probabilities stored in memory
	unsigned int n_storage_events;			// Number of storage having been done
	unsigned char curr_dipcodes [64];		// List of diplotypes in a given segment (buffer style variable)
	unsigned char curr_hapcodes [16];		// List of diplotypes in a given segment (buffer style variable)
	bool double_precision;


	// VARIANT / HAPLOTYPE / DIPLOTYPE DATA
	std::vector < unsigned char > Variants;		// 0.5 byte per variant
	std::vector < unsigned char > Ambiguous;		// 1 byte per ambiguous variant
	std::vector < unsigned long > Diplotypes;	// 8 bytes per segment
	std::vector < unsigned short > Lengths;		// 2 bytes per segment

	//PHASE PROBS
	std::vector < bool > ProbMask;
	std::vector < float > ProbStored;
	std::vector < float > ProbMissing;

	//METHODS
	genotype(unsigned int);
	~genotype();
	void free();
	void make(std::vector < unsigned char > &, std::vector < float > &);
	void make(std::vector < unsigned char > &);
	void build();
	void sample(std::vector < double > &, std::vector < float > &);
	void sampleForward(std::vector < double > &, std::vector < float > &);
	void sampleBackward(std::vector < double > &, std::vector < float > &);
	void solve();
	void mapMerges(std::vector < double > &, double , std::vector < bool > &);
	void performMerges(std::vector < double > &, std::vector < bool > &);
	void store(std::vector < double > &, std::vector < float > &);
	void scaffoldTrio(genotype *, genotype *, std::vector < unsigned int > &);
	void scaffoldDuoFather(genotype *, std::vector < unsigned int > &);
	void scaffoldDuoMother(genotype *, std::vector < unsigned int > &);

	//INLINES
	unsigned int countDiplotypes(unsigned long);
	void makeDiplotypes(unsigned long);
	unsigned int countTransitions();
	bool isOrdered(unsigned long _dip);
};

inline
bool genotype::isOrdered(unsigned long _dip) {
    std::fill(std::begin(curr_hapcodes), std::begin(curr_hapcodes)+16, 0);
	for (unsigned int d = 0, i = 0 ; d < 64 ; ++d) {
		if (DIP_GET(_dip, d)) {
			unsigned char hap0 = DIP_HAP0(d);
			unsigned char hap1 = DIP_HAP1(d);
			curr_hapcodes[hap0] = 1;
			curr_hapcodes[HAP_NUMBER + hap1] = 1;
		}
	}
	for (int h = 0 ; h < HAP_NUMBER ; h++) {
		if (curr_hapcodes[h] && curr_hapcodes[HAP_NUMBER+h]) {
			return false;
		}
	}
	return true;
}

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
