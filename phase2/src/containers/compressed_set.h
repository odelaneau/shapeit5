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
#ifndef _COMPRESSED_SET_H
#define _COMPRESSED_SET_H

#include <utils/otools.h>

struct cprob {
	unsigned int kst : 16;
	unsigned int lpb : 8;
	unsigned int rpb : 8;

	cprob(unsigned short _kst, unsigned short _lpb, unsigned short _rpb) {
		kst = _kst;
		lpb = _lpb;
		rpb = _rpb;
	}
};

class compressed_set {
public:

	unsigned int n_haplotypes;					//#samples

	vector < vector < cprob > > Phap_states;
	vector < vector < unsigned int > > Phap_indexes;

	//
	compressed_set();
	~compressed_set();
	void clear();
	void allocate(unsigned int );
};

#endif
