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

#include <containers/ibd2_tracks.h>

using namespace std;

ibd2_tracks::ibd2_tracks () {
	clear();
}

ibd2_tracks::~ibd2_tracks() {
	clear();
}

void ibd2_tracks::clear() {
	IBD2.clear();
}

void ibd2_tracks::initialize(int n_ind, variant_map & V) {
	IBD2 = vector < vector < track > > (n_ind);
	vec_cm = vector < float > (V.size(), 0.0f);
	for (int l = 0 ; l < vec_cm.size() ; l ++) vec_cm[l] = V.vec_pos[l]->cm;
}

void ibd2_tracks::pushIBD2(int ind, vector < track > & T) {
	expand(T);
	for (int t = 0 ; t < T.size() ; t ++)
		IBD2[min(ind, T[t].ind)].emplace_back(max(ind, T[t].ind), T[t].from, T[t].to);
}

void ibd2_tracks::expand(vector < track > & IBD) {
	unsigned n_merge = 0;
	for (int t = 0 ; t < IBD.size() ; t ++) {

		//Expand left by 4cm
		float leftcm = vec_cm[IBD[t].from];
		while (IBD[t].from > 1 && (leftcm - vec_cm[IBD[t].from]) < 4.0f) IBD[t].from --;

		//Expand right by 4cm
		float rightcm = vec_cm[IBD[t].to];
		while (IBD[t].to < (vec_cm.size()-1) && (vec_cm[IBD[t].to] - rightcm) < 4.0f) IBD[t].to ++;

	}
}

int ibd2_tracks::collapse(vector < track > & IBD) {
	unsigned n_merge = 0;
	if (IBD.size() > 1) {
		for(vector < track > :: iterator it = IBD.begin() + 1 ; it != IBD.end() ; ) {
			if (it->overlap(*(it-1))) {
				n_merge += (it-1)->merge(*it);
				it = IBD.erase(it);
			} else ++it;
		}
	}
	return n_merge;
}

void ibd2_tracks::collapse() {
	tac.clock();
	unsigned int n_inds1 = 0, n_tracks1 = 0, n_merged1 = 0, n_inds2 = 0, n_tracks2 = 0, n_merged2 = 0;
	for (int i = 0 ; i < IBD2.size() ; i ++) {
		sort(IBD2[i].begin(), IBD2[i].end());
		n_merged2 += collapse(IBD2[i]);
		n_tracks2 += IBD2[i].size();
		n_inds2 += (IBD2[i].size() > 0);
	}
	vrb.bullet("IBD2 tracks [#inds=" + stb.str(n_inds2) + " / #tracks=" + stb.str(n_tracks2) + " / #merged = " + stb.str(n_merged2) + "]");
}

bool ibd2_tracks::noIBD2(int hap0, int hap1, int locus) {
	int src_ind = min(hap0/2, hap1/2);
	int tar_ind = max(hap0/2, hap1/2);
	if (src_ind == tar_ind) return false;
	for (int i = 0 ; i < IBD2[src_ind].size() && IBD2[src_ind][i].ind <= tar_ind ; i ++)
		if ((IBD2[src_ind][i].ind == tar_ind) && (IBD2[src_ind][i].from <= locus) && (locus <= IBD2[src_ind][i].to)) return false;
	return true;
}




