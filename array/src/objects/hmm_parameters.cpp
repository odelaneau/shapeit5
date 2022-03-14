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

#include <objects/hmm_parameters.h>

hmm_parameters::hmm_parameters() {
	ed = 0.0001f;
	ee = 0.9999f;
}

hmm_parameters::~hmm_parameters() {
	t.clear();
	nt.clear();
}

void hmm_parameters::initialise(variant_map & V, int _Neff, int _Nhap) {
	Neff = _Neff; Nhap = _Nhap;
	cm = vector < float > (V.size(), 0.0);
	for (int l = 0 ; l < V.size() ; l ++) cm[l] = V.vec_pos[l]->cm;
	t = vector < float > (V.size() - 1, 0.0);
	nt = vector < float > (V.size() - 1, 0.0);
	for (int l = 1 ; l < cm.size() ; l ++) {
		float dist_cm = cm[l] - cm[l-1];
		if (dist_cm <= 1e-7) dist_cm = 1e-7;
		t[l-1] = -1.0f * expm1f(-0.04 * Neff * dist_cm / Nhap);
		nt[l-1] = 1-t[l-1];
	}
	int count_rare = 0;
	rare_allele = vector < char > (V.size(), -1);
	for (int l = 0 ; l < V.size() ; l ++) if (V.vec_pos[l]->getMAF() <= RARE_VARIANT_FREQ) {
		rare_allele[l] = (V.vec_pos[l]->getAF() > 0.5f);
		count_rare ++;
	}
	vrb.bullet("HMM parameters [Ne=" + stb.str(Neff) + " / Error=" + stb.str(ed) + " / #rare=" + stb.str(count_rare) + "]");
}

float hmm_parameters::getForwardTransProb(int prev_idx, int curr_idx) {
	assert(curr_idx>prev_idx);
	if (curr_idx == (prev_idx + 1)) return t[prev_idx];
	else {
		float dist_cm = cm[curr_idx] - cm[prev_idx];
		if (dist_cm <= 1e-7) dist_cm = 1e-7;
		return -1.0f * expm1f(-0.04 * Neff * dist_cm / Nhap);
	}
}

float hmm_parameters::getBackwardTransProb(int prev_idx, int curr_idx) {
	assert(curr_idx<prev_idx);
	if (curr_idx == (prev_idx - 1)) return t[curr_idx];
	else {
		float dist_cm = cm[prev_idx] - cm[curr_idx];
		if (dist_cm <= 1e-7) dist_cm = 1e-7;
		return -1.0f * expm1f(-0.04 * Neff * dist_cm / Nhap);
	}
}


