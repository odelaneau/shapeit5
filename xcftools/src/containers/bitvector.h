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

#ifndef _BITVECTOR_H
#define _BITVECTOR_H

#include <utils/otools.h>

class bitvector {
public:
	uint64_t n_bytes, n_elements;
	char * bytes;

	bitvector();
	bitvector(uint32_t size);
	~bitvector();

	void allocate(uint32_t size);
	void set(uint32_t idx, bool bit);
	void set(bool bit);
	bool get(uint32_t idx);
};

inline
void bitvector::set(uint32_t idx, bool value) {
	uint32_t idx_byt = idx / 8;
	uint32_t idx_bit = idx % 8;
	char mask = ~(1 << (7 - idx_bit));
	this->bytes[idx_byt] &= mask;
	this->bytes[idx_byt] |= (value << (7 - idx_bit));
}

inline
void bitvector::set(bool value) {
	if (value) memset(bytes, 0, n_bytes);
	else memset(bytes, UCHAR_MAX, n_bytes);
}

inline
bool bitvector::get(uint32_t idx) {
	uint32_t idx_byt = idx / 8;
	uint32_t idx_bit = idx % 8;
	return (this->bytes[idx_byt] >> (7 - (idx_bit%8))) & 1;
}

#endif
