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

#include <containers/bitvector.h>

using namespace std;

bitvector::bitvector() {
	n_bytes = n_elements = 0;
	bytes = NULL;
}

bitvector::bitvector(uint32_t size) {
	n_bytes = DIVU(size, 8);
	n_elements = size;
	bytes = (char*)malloc(n_bytes * sizeof(char));
	memset(bytes, 0, n_bytes);
}

bitvector::~bitvector() {
	n_bytes = n_elements = 0;
	if (bytes != NULL) free(bytes);
	bytes = NULL;
}

void bitvector::allocate(uint32_t size) {
	n_bytes = DIVU(size, 8);
	n_elements = size;
	bytes = (char*)malloc(n_bytes * sizeof(char));
	memset(bytes, 0, n_bytes);
}
