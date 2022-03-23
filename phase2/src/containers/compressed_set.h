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

#define FORCE_CAST(var, type) *(type*)&var

struct prob {
	unsigned int kst : 16;
	unsigned int lpb : 8;
	unsigned int rpb : 8;
};

class cprobs {
public:
	vector < unsigned int > values;

	cprobs(unsigned int _size, int _idx0, int _idx1) {
		values.reserve(_size+2);
		values.push_back(_idx0);
		values.push_back(_idx1);
	}

	~cprobs() {
		values.clear();
	}

	void push(unsigned short _kst, unsigned short _lpb, unsigned short _rpb) {
		prob p; p.kst = _kst; p.lpb = _lpb; p.rpb = _rpb;
		values.push_back(FORCE_CAST(p, unsigned int));
	}

	unsigned int size() {
		return values.size() - 2;
	}

	bool operator < (const cprobs & cp) const {
		return ((values[0]<cp.values[0]) || ((values[0]==cp.values[0])&&(values[1]<cp.values[1])));
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
