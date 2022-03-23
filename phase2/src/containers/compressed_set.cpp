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
	Pstates.clear();
}

void compressed_set::transpose() {
	tac.clock();
	for (unsigned long int e = 0 ; e < Pstates.size() ; e ++) swap(Pstates[e].idx0, Pstates[e].idx1);
	sort(Pstates.begin(), Pstates.end());
	vrb.bullet("Transpose compressed probabilities (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

unsigned long int compressed_set::sizeBytes() {
	unsigned long int nbytes = 0;
	for (unsigned long int e = 0 ; e < Pstates.size() ; e ++) {
		nbytes += 20;		//2 x unsigned int
		nbytes += Pstates[e].probs.size() * 4;
	}
	return nbytes;
}

void compressed_set::mapping(unsigned int n_scaffold_variants) {
	tac.clock();
	Pmapping = vector < long int > (n_scaffold_variants+1, -1);
	for (unsigned long int e = 0 ; e < Pstates.size() ; e ++) {
		if (Pmapping[Pstates[e].idx0] < 0) Pmapping[Pstates[e].idx0] = e;
	}
/*
	for (unsigned long int e = 0 ; e < Pstates.size() && e < 200; e ++) {
		cout << e << " " << Pstates[e].idx0 << " " << Pstates[e].idx1 << endl;
	}
*/
	vrb.bullet("Remap compressed probabilities (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
