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

#ifndef _BITMATRIX_H
#define _BITMATRIX_H

#include <utils/otools.h>

inline static uint32_t abracadabra(const uint32_t &i1, const uint32_t &i2) {
	return static_cast<unsigned int>((static_cast<uint64_t>(i1) * static_cast<uint64_t>(i2)) >> 32);
}

class bitmatrix	{
public:
	uint64_t n_bytes, n_cols, n_rows, startAddr;
	unsigned char * bytes;

	bitmatrix();
	~bitmatrix();

	void allocate(uint32_t nrow, uint32_t ncol);
	void allocateFast(uint32_t nrow, uint32_t ncol);

	void reallocateFull(uint32_t nrow, uint32_t ncol);
	void reallocateFast(uint32_t nrow, uint32_t ncol);
	void reallocate(uint32_t nrow, uint32_t ncol);


	void subset(bitmatrix & BM, std::vector < uint32_t > & rows);
	void getMatchHetCount(uint32_t i0, uint32_t i1, int32_t & c1, int32_t & m1);
	void set(uint32_t row, uint32_t col, unsigned char bit);
	unsigned char get(uint32_t row, uint32_t col);
	unsigned char getByte(uint32_t row, uint32_t col);
	void transpose(bitmatrix & BM, uint32_t _max_row, uint32_t _max_col);
	void transpose(bitmatrix & BM);
};

inline
void bitmatrix::set(uint32_t row, uint32_t col, unsigned char bit) {
	uint32_t bitcol = col % 8;
	uint64_t targetAddr = ((unsigned long)row) * (n_cols/8) + col/8;
	unsigned char mask = ~(1 << (7 - bitcol));
	this->bytes[targetAddr] &= mask;
	this->bytes[targetAddr] |= (bit << (7 - bitcol));
}

inline
unsigned char bitmatrix::get(uint32_t row, uint32_t col) {
	uint64_t targetAddr = ((unsigned long)row) * (n_cols>>3) +  (col>>3);
	return (this->bytes[targetAddr] >> (7 - (col%8))) & 1;
}

inline
unsigned char bitmatrix::getByte(uint32_t row, uint32_t col) {
	uint64_t targetAddr = ((unsigned long)row) * (n_cols>>3) +  (col>>3);
	return bytes[targetAddr];
}


#endif
