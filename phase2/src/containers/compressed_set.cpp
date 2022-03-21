////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <containers/compressed_set.h>

compressed_set::compressed_set() {
	clear();
}

compressed_set::~compressed_set() {
	clear();
}

void compressed_set::clear() {
	Phap_states.clear();
	Phap_indexes.clear();
}

void compressed_set::allocate(unsigned int _n_haplotypes) {
	tac.clock();
	n_haplotypes = _n_haplotypes;
	Phap_states = vector < vector < cprob > > (n_haplotypes);
	Phap_indexes = vector < vector < unsigned int > > (n_haplotypes);
	vrb.bullet("Compressed set allocation [#haplotypes=" + stb.str(n_haplotypes) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
