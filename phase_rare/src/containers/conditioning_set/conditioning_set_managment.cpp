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

#include <containers/conditioning_set/conditioning_set_header.h>

using namespace std;

conditioning_set::conditioning_set() {
	depth_common = 0;
	depth_rare = 0;
}

conditioning_set::~conditioning_set() {
	depth_common = 0;
	depth_rare = 0;
	sites_pbwt_evaluation.clear();
	sites_pbwt_selection.clear();
	sites_pbwt_grouping.clear();
	for (int h = 0 ; h < indexes_pbwt_neighbour.size() ; h ++) {
		indexes_pbwt_neighbour[h].clear();
		indexes_pbwt_neighbour[h].shrink_to_fit();
	}
	indexes_pbwt_neighbour.clear();
	indexes_pbwt_neighbour.shrink_to_fit();
}

void conditioning_set::initialize(variant_map & V, float _modulo_selection, float _mdr, int _depth_common, int _depth_rare, int _mac) {
	tac.clock();

	//SETTING PARAMETERS
	depth_common = _depth_common;
	depth_rare = _depth_rare;

	//MAPPING EVAL+GRP
	int n_evaluated = 0;
	sites_pbwt_evaluation = vector < bool > (V.sizeScaffold(), false);
	sites_pbwt_grouping = vector < int > (V.sizeScaffold(), -1);
	for (int l = 0 ; l < V.sizeScaffold() ; l ++) {
		sites_pbwt_evaluation[l] = (V.vec_scaffold[l]->getMAC() >= _mac && V.vec_scaffold[l]->getMDR() <= _mdr);
		sites_pbwt_grouping[l] = (int)round(V.vec_scaffold[l]->cm / _modulo_selection);
		n_evaluated += sites_pbwt_evaluation[l];
	}
	for (int l = 0, src = -1, tar = -1 ; l < V.sizeScaffold() ; l ++) {
		if (src == sites_pbwt_grouping[l]) sites_pbwt_grouping[l] = tar;
		else { src = sites_pbwt_grouping[l]; sites_pbwt_grouping[l] = ++tar; }
	}
	sites_pbwt_ngroups = sites_pbwt_grouping.back() + 1;

	//ALLOCATE
	indexes_pbwt_neighbour = vector < vector < unsigned int > > (n_haplotypes);
	IBD2 = vector < vector < unsigned int > > (n_haplotypes / 2);
	shuffledI = 0;
	shuffledO = vector < unsigned int > (n_haplotypes);
	iota(shuffledO.begin(), shuffledO.end(), 0);
	vrb.bullet("PBWT initialization [#eval=" + stb.str(n_evaluated) + " / #select=" + stb.str(sites_pbwt_grouping.back() + 1) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
