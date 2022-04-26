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
#include <containers/state_set.h>

state_set::state_set() {
	clear();
}

state_set::~state_set() {
	clear();
}

void state_set::clear() {
	Pstates1.clear();
	Pstates2.clear();
	Pmapping.clear();
	Pstates1.shrink_to_fit();
	Pstates2.shrink_to_fit();
	Pmapping.shrink_to_fit();
}

void state_set::transpose1() {
	tac.clock();
	for (unsigned long int e = 0 ; e < Pstates1.size() ; e ++) Pstates1[e].swap();
	sort(Pstates1.begin(), Pstates1.end());
	vrb.bullet("Transpose1 compressed probabilities (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void state_set::transpose2() {
	tac.clock();
	for (unsigned long int e = 0 ; e < Pstates2.size() ; e ++) Pstates2[e].swap();
	sort(Pstates2.begin(), Pstates2.end());
	vrb.bullet("Transpose2 compressed probabilities (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void state_set::mapping1(unsigned int n_scaffold_variants) {
	tac.clock();
	Pmapping = vector < long int > (n_scaffold_variants+1, -1);
	for (unsigned long int e = 0 ; e < Pstates1.size() ; e ++) {
		if (Pmapping[Pstates1[e].id0] < 0) Pmapping[Pstates1[e].id0] = e;
	}
	vrb.bullet("Map1 compressed probabilities (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void state_set::mapping2(unsigned int n_scaffold_variants) {
	tac.clock();
	Pmapping = vector < long int > (n_scaffold_variants+1, -1);
	for (unsigned long int e = 0 ; e < Pstates2.size() ; e ++) {
		if (Pmapping[Pstates2[e].id0] < 0) Pmapping[Pstates2[e].id0] = e;
	}
	vrb.bullet("Map2 compressed probabilities (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
