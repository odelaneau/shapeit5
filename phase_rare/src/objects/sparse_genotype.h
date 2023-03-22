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

#ifndef _SPARSE_GENOTYPE_H
#define _SPARSE_GENOTYPE_H

#include <utils/otools.h>

#define SETBIT(n,i)	(n)|=(1U<<(i));
#define CLRBIT(n,i)	(n)&=~(1U<<i);
#define GETBIT(n,i)	(((n)>>(i))&1U);


class sparse_genotype {
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

	sparse_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
		prob = -1.0f;
	}

	sparse_genotype(unsigned int value) {
		set(value);
		prob = pha?1.0f:-1.0f;
	}

	sparse_genotype(unsigned int _idx, bool _het, bool _mis, bool _al0, bool _al1, bool _pha) {
		idx = _idx; het = _het; mis = _mis; al0 = _al0; al1 = _al1;

		pha = _pha || (!het && !mis);

		if (pha) prob = 1.0f;
		else {
			prob = -1.0f;
			if (al0 != al1) {
				if (rng.flipCoin()) { al0 = 0; al1 = 1; }
				else { al0 = 1; al1 = 0; }
			}
		}
	}

	~sparse_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
		prob = -1.0f;
	}

	bool operator < (const sparse_genotype & rg) const {
		return idx < rg.idx;
	}

	int get(bool majorA) {
		if (mis) return -1;
		if (het) return 1;
		return 2-majorA;
	}

	unsigned int get() {
		unsigned int value = (idx << 5);
		if (het) SETBIT(value, 4);
		if (mis) SETBIT(value, 3);
		if (al0) SETBIT(value, 2);
		if (al1) SETBIT(value, 1);
		if (pha) SETBIT(value, 0);
		return value;
	}

	void set(unsigned int value) {
		idx = (value >> 5);
		het = GETBIT(value, 4);
		mis = GETBIT(value, 3);
		al0 = GETBIT(value, 2);
		al1 = GETBIT(value, 1);
		pha = GETBIT(value, 0);
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
			std::vector < double > gprobs = std::vector < double > (4, 0.0f);
			float p01 = std::max(prb0, std::numeric_limits<float>::min());
			float p00 = std::max(1.0f - prb0, std::numeric_limits<float>::min());
			float p11 = std::max(prb1, std::numeric_limits<float>::min());
			float p10 = std::max(1.0f - prb1, std::numeric_limits<float>::min());

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

			int maxg = alg.imax(gprobs);
			switch (maxg) {
			case 0:	al0 = 0; al1 = 0; break;
			case 1:	al0 = 0; al1 = 1; break;
			case 2:	al0 = 1; al1 = 0; break;
			case 3:	al0 = 1; al1 = 1; break;
			}
			prob = gprobs[maxg] / (gprobs[0] + gprobs[1] + gprobs[2] + gprobs[3]);
		}
	}

	std::string str() {
		return std::string("[" + stb.str(idx) + "/" + stb.str(het) + "/" + stb.str(mis) + "/" + stb.str(al0) + "/" + stb.str(al1) + "/" + stb.str(pha) + "]");
	}
};

#endif
