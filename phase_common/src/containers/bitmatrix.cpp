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

#include <containers/bitmatrix.h>

using namespace std;

static unsigned char nbit_set[256] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };

bitmatrix::bitmatrix() {
	n_rows = 0;
	n_cols = 0;
	n_bytes = 0;
	bytes = NULL;
}

bitmatrix::~bitmatrix() {
	n_bytes = 0;
	if (bytes != NULL) free(bytes);
}

int bitmatrix::subset(bitmatrix & BM, vector < unsigned int > rows, unsigned int col_from, unsigned int col_to) {
	n_rows = rows.size() + ((rows.size()%8)?(8-(rows.size()%8)):0);
	unsigned long row_start = col_from/8;
	unsigned long row_end = col_to/8;
	unsigned long n_bytes_per_row = row_end - row_start + 1;
	n_cols = n_bytes_per_row * 8;
	n_bytes = n_bytes_per_row * n_rows;
	bytes = (unsigned char*)malloc(n_bytes*sizeof(unsigned char));
	unsigned long offset_addr = 0;
	for (int r = 0 ; r < rows.size() ; r ++) {
		row_start = ((unsigned long)rows[r]) * (BM.n_cols/8) + col_from/8;
		row_end = ((unsigned long)rows[r]) * (BM.n_cols/8) + col_to/8;
		memcpy(&bytes[offset_addr], &BM.bytes[row_start], n_bytes_per_row);
		offset_addr += n_bytes_per_row;
	}
	return col_from % 8;
}

/*
void bitmatrix::getMatchHetCount(unsigned int i0, unsigned int i1, unsigned int start, unsigned int stop, int & c1, int & m1) {
	c1=m1=0;
	unsigned long n_bytes_per_row = stop/8 - start/8 + 1;
	unsigned long offset_i0_h0 = (unsigned long)(2*i0+0)*(n_cols/8) + start/8;
	unsigned long offset_i0_h1 = (unsigned long)(2*i0+1)*(n_cols/8) + start/8;
	unsigned long offset_i1_h0 = (unsigned long)(2*i1+0)*(n_cols/8) + start/8;
	unsigned long offset_i1_h1 = (unsigned long)(2*i1+1)*(n_cols/8) + start/8;
	for (unsigned long b = 0 ; b < n_bytes_per_row ; b ++) {
		//unsigned char i0_g0 = (bytes[offset_i0_h0+b]|bytes[offset_i0_h1+b]);
		unsigned char i0_g1 = (bytes[offset_i0_h0+b]^bytes[offset_i0_h1+b]);
		//unsigned char i0_g2 = (bytes[offset_i0_h0+b]&bytes[offset_i0_h1+b]);
		//unsigned char i1_g0 = (bytes[offset_i1_h0+b]|bytes[offset_i1_h1+b]);
		unsigned char i1_g1 = (bytes[offset_i1_h0+b]^bytes[offset_i1_h1+b]);
		//unsigned char i1_g2 = (bytes[offset_i1_h0+b]&bytes[offset_i1_h1+b]);
		//m0 += distance_lookup[i0_g0 ^ i1_g0];
		m1 += nbit_set[i0_g1 ^ i1_g1];
		//m2 += distance_lookup[i0_g2 ^ i1_g2];
		c1 += nbit_set[i0_g1 | i1_g1];
	}
}

void bitmatrix::getMatchHetCount_seq(unsigned int i0, unsigned int i1, unsigned int start, unsigned int stop, int & _ibd2_start, int & _ibd2_stop) {
	unsigned long n_bytes_per_row = stop/8 - start/8 + 1;
	unsigned long offset_i0_h0 = (unsigned long)(2*i0+0)*(n_cols/8) + start/8;
	unsigned long offset_i0_h1 = (unsigned long)(2*i0+1)*(n_cols/8) + start/8;
	unsigned long offset_i1_h0 = (unsigned long)(2*i1+0)*(n_cols/8) + start/8;
	unsigned long offset_i1_h1 = (unsigned long)(2*i1+1)*(n_cols/8) + start/8;

	int maxZeroCount = 0;
	int maxStartIdx = -1;
	int maxStopIdx = -1;
	int currZeroCount = 0;
	int currStartIdx = -1;
	int currStopIdx = -1;

	for (unsigned long b = 0 ; b < n_bytes_per_row ; b ++) {
		unsigned char i0_g1 = (bytes[offset_i0_h0+b]^bytes[offset_i0_h1+b]);
		unsigned char i0_g2 = (bytes[offset_i0_h0+b]&bytes[offset_i0_h1+b]);
		unsigned char i1_g1 = (bytes[offset_i1_h0+b]^bytes[offset_i1_h1+b]);
		unsigned char i1_g2 = (bytes[offset_i1_h0+b]&bytes[offset_i1_h1+b]);
		bool diff1 = ((i0_g1 ^ i1_g1) != 0);
		bool diff2 = ((i0_g2 ^ i1_g2) != 0);

		if (diff1 || diff2) currZeroCount = 0;
		else {
			if (currZeroCount == 0) {
				currStartIdx = b;
				currStopIdx = b;
			} else currStopIdx++;
			currZeroCount++;
			if (currZeroCount > maxZeroCount) {
				maxZeroCount = currZeroCount;
				maxStartIdx = currStartIdx;
				maxStopIdx = currStopIdx;
			}
		}
	}

	if (maxZeroCount > 0) {
		_ibd2_start = max(maxStartIdx*8 + start, start);
		_ibd2_stop = min(maxStopIdx*8 + start, stop);
	} else {
		_ibd2_start = -1;
		_ibd2_stop = -1;
	}
}

void bitmatrix::getMatchHets(unsigned int i0, unsigned int i1, unsigned int start, unsigned int stop, float & fine_proportion, float & approx_proportion, float & transitionS2S, float & transitionS2D, float & transitionD2S, float & transitionD2D) {
	unsigned long n_bytes_per_row = stop/8 - start/8 + 1;
	unsigned long offset_i0_h0 = (unsigned long)(2*i0+0)*(n_cols/8) + start/8;
	unsigned long offset_i0_h1 = (unsigned long)(2*i0+1)*(n_cols/8) + start/8;
	unsigned long offset_i1_h0 = (unsigned long)(2*i1+0)*(n_cols/8) + start/8;
	unsigned long offset_i1_h1 = (unsigned long)(2*i1+1)*(n_cols/8) + start/8;

	bool lastSame = true;
	unsigned int interHet = 0, unionHet = 0;
	unsigned int interBloc = 0, unionBloc = 0;
	unsigned int transS2S = 0, transD2S = 0, transS2D = 0, transD2D = 0;
	for (unsigned long b = 0 ; b < n_bytes_per_row ; b ++) {
		unsigned char i0_g1 = (bytes[offset_i0_h0+b]^bytes[offset_i0_h1+b]);
		unsigned char i1_g1 = (bytes[offset_i1_h0+b]^bytes[offset_i1_h1+b]);

		unsigned int iHet = nbit_set[i0_g1 ^ i1_g1];	// iHet != 0 when non matching hets in the block
		unsigned int uHet = nbit_set[i0_g1 | i1_g1];	// uHet != 0 when hets are in the block

		interHet += iHet;
		unionHet += uHet;

		interBloc += (iHet>0);
		unionBloc += (uHet>0);

		if (uHet>0) {
			transS2S += (lastSame && !iHet);
			transS2D += (lastSame && iHet);
			transD2S += (!lastSame && !iHet);
			transD2D += (!lastSame && iHet);
			lastSame = !iHet;
		}
	}

	fine_proportion = (unionHet - interHet) * 1.0f / unionHet;
	approx_proportion = (unionBloc - interBloc) * 1.0f / unionBloc;

	float sum = transS2S + transD2S + transS2D + transD2D;
	transitionS2S = transS2S * 1.0f / sum;
	transitionS2D = transS2D * 1.0f / sum;
	transitionD2S = transD2S * 1.0f / sum;
	transitionD2D = transD2D * 1.0f / sum;
}
*/

float bitmatrix::getMatchHets(unsigned int i0, unsigned int i1, unsigned int start, unsigned int stop) {
	unsigned long n_bytes_per_row = stop/8 - start/8 + 1;
	unsigned long offset_i0_h0 = (unsigned long)(2*i0+0)*(n_cols/8) + start/8;
	unsigned long offset_i0_h1 = (unsigned long)(2*i0+1)*(n_cols/8) + start/8;
	unsigned long offset_i1_h0 = (unsigned long)(2*i1+0)*(n_cols/8) + start/8;
	unsigned long offset_i1_h1 = (unsigned long)(2*i1+1)*(n_cols/8) + start/8;

	unsigned int interHet = 0, unionHet = 0;
	for (unsigned long b = 0 ; b < n_bytes_per_row ; b ++) {
		unsigned char i0_g1 = (bytes[offset_i0_h0+b]^bytes[offset_i0_h1+b]);
		unsigned char i1_g1 = (bytes[offset_i1_h0+b]^bytes[offset_i1_h1+b]);
		interHet += nbit_set[i0_g1 ^ i1_g1];	// iHet != 0 when non matching hets in the block
		unionHet += nbit_set[i0_g1 | i1_g1];	// uHet != 0 when hets are in the block
	}
	return (unionHet)? ((unionHet - interHet) * 1.0f / unionHet) : 0.0f;
}


void bitmatrix::allocate(unsigned int nrow, unsigned int ncol) {
	n_rows = nrow + ((nrow%8)?(8-(nrow%8)):0);
	n_cols = ncol + ((ncol%8)?(8-(ncol%8)):0);
	n_bytes = (n_cols/8) * (unsigned long)n_rows;
	bytes = (unsigned char*)malloc(n_bytes*sizeof(unsigned char));
	memset(bytes, 0, n_bytes);
}

void bitmatrix::allocateFast(unsigned int nrow, unsigned int ncol) {
	n_rows = nrow + ((nrow%8)?(8-(nrow%8)):0);
	n_cols = ncol + ((ncol%8)?(8-(ncol%8)):0);
	n_bytes = (n_cols/8) * (unsigned long)n_rows;
	bytes = (unsigned char*)malloc(n_bytes*sizeof(unsigned char));
}


/*
 * This algorithm for transposing bit matrices is adapted from the code of Timur Kristóf
 * Timur Kristóf: https://github.com/venemo
 * Original version of the code (MIT license): https://github.com/Venemo/fecmagic/blob/master/src/binarymatrix.h
 * Of note, function abracadabra is the same than getMultiplyUpperPart function in the original code from Timur Kristóf.
 */
void bitmatrix::transpose(bitmatrix & BM, unsigned int _max_row, unsigned int _max_col) {
	unsigned int max_row = _max_row + ((_max_row%8)?(8-(_max_row%8)):0);
	unsigned int max_col = _max_col + ((_max_col%8)?(8-(_max_col%8)):0);
	unsigned long targetAddr, sourceAddr;
	union { unsigned int x[2]; unsigned char b[8]; } m4x8d;
	for (unsigned int row = 0; row < max_row; row += 8) {
		for (unsigned int col = 0; col < max_col; col += 8) {
			for (unsigned int i = 0; i < 8; i++) {
				sourceAddr = (row+i) * ((unsigned long)(n_cols/8)) + col/8;
				m4x8d.b[7 - i] = this->bytes[sourceAddr];
			}
			for (unsigned int i = 0; i < 7; i++) {
				targetAddr = ((col+i) * ((unsigned long)(n_rows/8)) + (row) / 8);
				BM.bytes[targetAddr]  = static_cast<unsigned char>(abracadabra(m4x8d.x[1] & (0x80808080 >> i), (0x02040810 << i)) & 0x0f) << 4;
                BM.bytes[targetAddr] |= static_cast<unsigned char>(abracadabra(m4x8d.x[0] & (0x80808080 >> i), (0x02040810 << i)) & 0x0f) << 0;
			}
			targetAddr = ((col+7) * ((unsigned long)(n_rows/8)) + (row) / 8);
			BM.bytes[targetAddr]  = static_cast<unsigned char>(abracadabra((m4x8d.x[1] << 7) & (0x80808080 >> 0), (0x02040810 << 0)) & 0x0f) << 4;
            BM.bytes[targetAddr] |= static_cast<unsigned char>(abracadabra((m4x8d.x[0] << 7) & (0x80808080 >> 0), (0x02040810 << 0)) & 0x0f) << 0;
		}
	}
}

void bitmatrix::transpose(bitmatrix & BM) {
	transpose(BM, n_rows, n_cols);
}
