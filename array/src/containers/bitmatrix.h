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
#ifndef _BITMATRIX_H
#define _BITMATRIX_H

#include <utils/otools.h>

inline static unsigned int abracadabra(const unsigned int &i1, const unsigned int &i2) {
	return static_cast<unsigned int>((static_cast<unsigned long int>(i1) * static_cast<unsigned long int>(i2)) >> 32);
}

class bitmatrix	{
public:
	unsigned long int n_bytes, n_cols, n_rows, startAddr;
	unsigned char * bytes;

	bitmatrix();
	~bitmatrix();

	int subset(bitmatrix & BM, vector < unsigned int > rows, unsigned int col_from, unsigned int col_to);
	void getMatchHetCount(unsigned int i0, unsigned int i1, unsigned int start, unsigned int stop, int & c1, int & m1);
	void allocate(unsigned int nrow, unsigned int ncol);
	void allocateFast(unsigned int nrow, unsigned int ncol);
	void set(unsigned int row, unsigned int col, unsigned char bit);
	unsigned char get(unsigned int row, unsigned int col);
	void transpose(bitmatrix & BM, unsigned int _max_row, unsigned int _max_col);
	void transpose(bitmatrix & BM);
};

inline
void bitmatrix::set(unsigned int row, unsigned int col, unsigned char bit) {
	unsigned int bitcol = col % 8;
	unsigned long targetAddr = ((unsigned long)row) * (n_cols/8) + col/8;
	unsigned char mask = ~(1 << (7 - bitcol));
	this->bytes[targetAddr] &= mask;
	this->bytes[targetAddr] |= (bit << (7 - bitcol));
}

inline
unsigned char bitmatrix::get(unsigned int row, unsigned int col) {
	unsigned long targetAddr = ((unsigned long)row) * (n_cols>>3) +  (col>>3);
	return (this->bytes[targetAddr] >> (7 - (col%8))) & 1;
}

#endif
