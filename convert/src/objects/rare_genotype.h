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

#ifndef _RARE_GENOTYPE_H
#define _RARE_GENOTYPE_H

#include <utils/otools.h>

#define SETBIT(n,i)	(n)|=(1UL<<(i));
#define CLRBIT(n,i)	(n)&=~(1UL<<i);
#define GETBIT(n,i)	(((n)>>(i))&1U);


class rare_genotype {
public:

	unsigned int idx : 27;	//Individual index of the non Major/Major genotype
	unsigned int het : 1;	//Is it het?
	unsigned int mis : 1;	//Is it missing?	If none of the two, it is then Minor/Minor
	unsigned int al0 : 1;	//What's the allele on haplotype 0?
	unsigned int al1 : 1;	//What's the allele on haplotype 1?
	unsigned int pha : 1;	//Is the genotype phased?

	rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
	}

	rare_genotype(unsigned int _idx, bool _het, bool _mis, bool _al0, bool _al1, bool _pha) {
		idx = _idx; het = _het; mis = _mis; al0 = _al0; al1 = _al1; pha = _pha;
	}

	~rare_genotype() {
		idx = het = mis = al0 = al1 = pha = 0;
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
};

#endif
