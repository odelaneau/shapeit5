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

struct prob {
	unsigned int kst : 16;
	unsigned int lpb : 8;
	unsigned int rpb : 8;

	prob(unsigned short _kst, unsigned short _lpb, unsigned short _rpb) {
		kst = _kst;
		lpb = _lpb;
		rpb = _rpb;
	}
};

class cprobs {
public:
	int idx0;
	int idx1;
	vector < prob > probs;

	cprobs(int _idx0, int _idx1) {
		idx0 = _idx0;
		idx1 = _idx1;
	}

	~cprobs() {
		idx0 = -1; idx1 = -1;
		probs.clear();
	}

	void push(unsigned short _kst, unsigned short _lpb, unsigned short _rpb) {
		probs.emplace_back(_kst, _lpb, _rpb);
	}

	unsigned int size() {
		return probs.size();
	}

	bool operator < (const cprobs & cp) const {
		return ((idx0<cp.idx0) || ((idx0==cp.idx0)&&(idx1<cp.idx1)));
	}
};

class compressed_set {
public:
	vector < cprobs > Pstates;
	vector < long int > Pmapping;

	//
	compressed_set();
	~compressed_set();
	void clear();
	void transpose();
	void mapping(unsigned int);

	unsigned long int sizeBytes();
};

#endif
