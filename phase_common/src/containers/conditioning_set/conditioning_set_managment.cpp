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
	depth = 0;
	nthread = 0;
}

conditioning_set::~conditioning_set() {
	if (nthread > 1) {
		i_worker = 0;
		i_job = 0;
		pthread_mutex_destroy(&mutex_workers);
		id_workers.clear();
	}
	depth = 0;
	nthread = 0;
	sites_pbwt_mthreading.clear();
	sites_pbwt_evaluation.clear();
	sites_pbwt_selection.clear();
	sites_pbwt_grouping.clear();
	indexes_pbwt_neighbour.clear();
}


bool conditioning_set::split(variant_map & V, float min_length, int left_index, int right_index, vector < int > & output) {
	int chunkSize = right_index - left_index + 1;
	float chunkLength = V.vec_pos[right_index]->cm - V.vec_pos[left_index]->cm;

	if ((chunkSize > 2) && (chunkLength > min_length)) {
		vector <  int > left_output, right_output;
		bool ret1 = split(V, min_length, left_index, left_index + chunkSize/2 - 1, left_output);
		bool ret2 = split(V, min_length, left_index + chunkSize/2, right_index, right_output);

		if (ret1 && ret2) {
			output = vector < int >(left_output.size() + right_output.size());
			std::copy(left_output.begin(), left_output.end(), output.begin());
			std::copy(right_output.begin(), right_output.end(), output.begin() + left_output.size());
		} else {
			output.clear();
			output.push_back(left_index);
			output.push_back(right_index);
		}
		return true;
	} else return false;
}

void conditioning_set::initialize(variant_map & V, float _modulo_selection, float _modulo_multithreading, float _mdr, int _depth, int _mac, int _nthread) {
	tac.clock();

	//SETTING PARAMETERS
	depth = _depth;
	nthread = _nthread;
	if (nthread > 1) {
		i_worker = 0;
		i_job = 0;
		id_workers = vector < pthread_t > (nthread);
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//MAPPING EVAL+GRP
	int n_evaluated = 0;
	sites_pbwt_evaluation = vector < bool > (V.size(), false);
	sites_pbwt_mthreading = vector < int > (V.size(), -1);
	sites_pbwt_grouping = vector < int > (V.size(), -1);
	for (int l = 0 ; l < V.size() ; l ++) {
		sites_pbwt_evaluation[l] = (V.vec_pos[l]->getMAC() >= _mac && V.vec_pos[l]->getMDR() <= _mdr);
		sites_pbwt_grouping[l] = (int)round(V.vec_pos[l]->cm / _modulo_selection);
		n_evaluated += sites_pbwt_evaluation[l];
	}
	for (int l = 0, src = -1, tar = -1 ; l < V.size() ; l ++) {
		if (src == sites_pbwt_grouping[l]) sites_pbwt_grouping[l] = tar;
		else { src = sites_pbwt_grouping[l]; sites_pbwt_grouping[l] = ++tar; }
	}
	sites_pbwt_ngroups = sites_pbwt_grouping.back() + 1;

	//MAPPING MT+STOR
	vector < int > outputMT; outputMT.push_back(0); outputMT.push_back(V.size() - 1);
	split(V, _modulo_multithreading, 0, V.size() - 1, outputMT);
	for (int c = 0, ts = 0 ; c < outputMT.size() ; c += 2) {
		for (int l = outputMT[c] ; l <= outputMT[c+1] ; l++) sites_pbwt_mthreading[l] = c/2;
	}

	//MAPPING START OF MT
	starts_pbwt_mthreading = vector < int > (sites_pbwt_mthreading.back() + 1, -1);
	for (int l = 0 ; l < V.size() ; l ++) if (starts_pbwt_mthreading[sites_pbwt_mthreading[l]] < 0) starts_pbwt_mthreading[sites_pbwt_mthreading[l]] = l;
	for (int c = 0 ; c < starts_pbwt_mthreading.size() ; c ++) {
		float distance_cm = 0.0f;
		float storage_pos = V.vec_pos[starts_pbwt_mthreading[c]]->cm;
		while (starts_pbwt_mthreading[c] > 0 && distance_cm < 0.5f) {
			starts_pbwt_mthreading[c] --;
			distance_cm = storage_pos - V.vec_pos[starts_pbwt_mthreading[c]]->cm;
		}
	}

	//ALLOCATE
	Kbanned.initialize(n_ind);
	indexes_pbwt_neighbour = vector < int > ( (depth+1) * (sites_pbwt_grouping.back()+1) * n_ind * 2UL, -1);
	vrb.bullet("PBWT initialization [#eval=" + stb.str(n_evaluated) + " / #select=" + stb.str(sites_pbwt_grouping.back() + 1) + " / #chunk=" + stb.str(sites_pbwt_mthreading.back() + 1) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
