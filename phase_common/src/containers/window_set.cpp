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

#include <containers/window_set.h>

using namespace std;

window_set::window_set() {
	W.clear();
}

window_set::~window_set() {
	W.clear();
}

void window_set::clear() {
	W.clear();
}

int window_set::size() {
	return W.size();
}

bool window_set::split(double min_length_cm, int left_index, int right_index, vector < int > & idx_sta, vector < int > & idx_sto, vector < double > & ccm_sta, vector < double > & ccm_sto, vector < int > & output) {
	int number_of_segments = right_index-left_index+1;
	int number_of_variants = idx_sto[right_index] - idx_sta[left_index] + 1;
	double length_of_region = ccm_sto[right_index] - ccm_sta[left_index];

	//A phasing window must (i) span >=4 segments, (ii) contain >= 100 variants and (iii) span more than "min_length_cm" cM
	if (number_of_segments < 4 || number_of_variants < 100 || length_of_region < min_length_cm) return false;
	else {
		int split_point = rng.getInt(number_of_segments/2) + number_of_segments/4 + 1;
		vector <  int > left_output, right_output;
		bool ret1 = split(min_length_cm, left_index, left_index + split_point, idx_sta, idx_sto, ccm_sta, ccm_sto, left_output);
		bool ret2 = split(min_length_cm, left_index + split_point, right_index, idx_sta, idx_sto, ccm_sta, ccm_sto, right_output);

		if (ret1 && ret2) {
			//succesful split, so operate it
			output = vector < int >(left_output.size() + right_output.size());
			std::copy(left_output.begin(), left_output.end(), output.begin());
			std::copy(right_output.begin(), right_output.end(), output.begin() + left_output.size());
		} else {
			//unsuccesful split, so return current coordinates
			output.clear();
			output.push_back(left_index);
			output.push_back(right_index);
		}
		return true;
	}
}


int window_set::build (variant_map & V, genotype * g, float min_window_size) {

	//1. Mapping coordinates of each segment
	vector < unsigned int > loc_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > loc_siz = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > amb_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > amb_siz = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > mis_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > mis_siz = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > tra_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > tra_siz = vector < unsigned int >(g->n_segments, 0);
	vector < double > ccm_sta = vector < double >(g->n_segments, 0);
	vector < double > ccm_sto = vector < double >(g->n_segments, 0);
	vector < int > idx_sta = vector < int >(g->n_segments, 0);
	vector < int > idx_sto = vector < int >(g->n_segments, 0);

	unsigned int prev_dipcounts = 1, curr_dipcounts = 0;
	for (unsigned int s = 0, a = 0, t = 0, v = 0, m = 0 ; s < g->n_segments ; s ++) {
		//update a
		amb_idx[s] = a;
		for (unsigned int vrel = 0 ; vrel < g->Lengths[s] ; vrel ++) amb_siz[s] += VAR_GET_AMB(MOD2(v+vrel), g->Variants[DIV2(v+vrel)]);
		a += amb_siz[s];
		//update m
		mis_idx[s] = m;
		for (unsigned int vrel = 0 ; vrel < g->Lengths[s] ; vrel ++) mis_siz[s] += VAR_GET_MIS(MOD2(v+vrel), g->Variants[DIV2(v+vrel)]);
		m += mis_siz[s];
		//update v
		loc_idx[s] = v;
		loc_siz[s] = g->Lengths[s];
		v += loc_siz[s];
		//update idx
		idx_sta[s] = loc_idx[s];
		idx_sto[s] = loc_idx[s]+loc_siz[s]-1;
		//update ccm
		ccm_sta[s] = V.vec_pos[idx_sta[s]]->cm;
		ccm_sto[s] = V.vec_pos[idx_sto[s]]->cm;
		//update t
		tra_idx[s] = t;
		curr_dipcounts = g->countDiplotypes(g->Diplotypes[s]);
		tra_siz[s] = prev_dipcounts * curr_dipcounts;
		t += tra_siz[s];
		prev_dipcounts = curr_dipcounts;
	}

	//2. Reccursive split
	vector < int > output;
	output.push_back(0);
	output.push_back(g->n_segments-1);
	split(min_window_size, 0, g->n_segments-1, idx_sta, idx_sto, ccm_sta, ccm_sto, output);
	int n_windows = output.size()/2;

	//3. Update coordinates
	W = vector < window > (n_windows);
	for (unsigned int w = 0 ; w < n_windows ; w ++) {
		W[w].start_segment = output[2*w+0];
		W[w].stop_segment = output[2*w+1];
		W[w].start_ambiguous = amb_idx[W[w].start_segment];
		W[w].stop_ambiguous = amb_idx[W[w].stop_segment] + amb_siz[W[w].stop_segment] - 1;
		W[w].start_missing = mis_idx[W[w].start_segment];
		W[w].stop_missing = mis_idx[W[w].stop_segment] + mis_siz[W[w].stop_segment] - 1;
		W[w].start_locus = loc_idx[W[w].start_segment];
		W[w].stop_locus = loc_idx[W[w].stop_segment] + loc_siz[W[w].stop_segment] - 1;
		W[w].start_transition = tra_idx[W[w].start_segment] + tra_siz[W[w].start_segment];
		W[w].stop_transition = tra_idx[W[w].stop_segment] + tra_siz[W[w].stop_segment] - 1;
	}
	return n_windows;
}
