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

#ifndef _STATE_SET_H
#define _STATE_SET_H

#include <utils/otools.h>

class cstate {
public:
	uint32_t id0;
	uint32_t id1;
	uint32_t kst : 16;
	uint32_t lpb : 8;
	uint32_t rpb : 8;

	cstate(uint32_t _id0, uint32_t _id1, uint32_t _kst, unsigned char _lpb, unsigned char _rpb) {
		id0 = _id0; id1 = _id1; kst = _kst; lpb = _lpb; rpb = _rpb;
	}

	~cstate(){
	}

	bool operator < (const cstate & s) const {
		return ((id0<s.id0) || ((id0==s.id0)&&(id1<s.id1)));
	}

	void swap() {
		uint32_t tmp = id1; id1 = id0; id0 = tmp;
	}
};

class state_set {
public:
	std::vector < cstate > Pstates;
	std::vector < int64_t > Pmapping;

	//
	state_set();
	~state_set();
	void clear();
	void transpose();
	void mapping(uint32_t);
};

#endif
