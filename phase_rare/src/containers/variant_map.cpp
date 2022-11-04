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

#include <containers/variant_map.h>

using std::vector;
using std::pair;
using std::multimap;

variant_map::variant_map() {
}

variant_map::~variant_map() {
	for (int s = 0 ; s < vec_full.size() ; s++) delete vec_full[s];
	vec_full.clear();
	vec_common.clear();
	map_pos.clear();
}

unsigned int variant_map::sizeFull() {
	return vec_full.size();
}

unsigned int variant_map::sizeCommon() {
	return vec_common.size();
}

unsigned int variant_map::sizeRare() {
	return vec_rare.size();
}

unsigned int variant_map::sizeScaffold() {
	return vec_scaffold.size();
}

void variant_map::getCommonVariants(unsigned int vs0, unsigned int vs1, vector < unsigned int > & VC) {
	assert(vs1>vs0);
	VC.clear();
	unsigned int idx_vf0 = vec_scaffold[vs0]->idx_full;
	unsigned int idx_vf1 = vec_scaffold[vs1]->idx_full;
	for (int v = idx_vf0 ; v <= idx_vf1 ; v++)
		if (vec_full[v]->type == VARTYPE_COMM)
			VC.push_back(vec_full[v]->idx_common);
}

void variant_map::getRareVariants(unsigned int vs0, unsigned int vs1, vector < unsigned int > & VR) {
	assert(vs1>vs0);
	VR.clear();
	unsigned int idx_vf0 = vec_scaffold[vs0]->idx_full;
	unsigned int idx_vf1 = vec_scaffold[vs1]->idx_full;
	for (int v = idx_vf0 ; v <= idx_vf1 ; v++)
		if (vec_full[v]->type == VARTYPE_RARE)
			VR.push_back(vec_full[v]->idx_rare);
}

vector < variant * > variant_map::getByPos (int pos) {
	vector < variant * > vecS;
	pair < multimap < int , variant * >::iterator , multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}

void variant_map::push(variant * v) {
	if (v->type == VARTYPE_SCAF) {
		v->idx_scaffold = vec_scaffold.size();
		vec_scaffold.push_back(v);
	} else if (v->type == VARTYPE_COMM) {
		v->idx_common = vec_common.size();
		vec_common.push_back(v);
	} else {
		v->idx_rare = vec_rare.size();
		vec_rare.push_back(v);
	}
	v->idx_full = vec_full.size();
	vec_full.push_back(v);
	map_pos.insert(pair < int , variant * > (v->bp, v));
}

int variant_map::setCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM) {
	int cpt = 0;
	for (int l = 0 ; l < pos_cM.size() ; l ++) {
		vector  < variant * > vecS = getByPos(pos_bp[l]);
		for (int si = 0 ; si < vecS.size() ; si ++) {
			vecS[si]->cm = pos_cM[l];
			cpt++;
		}
	}
	return cpt;
}

int variant_map::interpolateCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM) {
	int n_interpolated = 0, i_locus = 0;
	double base, rate, dist;
	double mean_rate = (pos_cM.back() - pos_cM[0]) / (pos_bp.back() - pos_bp[0]);

	//Set up first positions to be mean rate
	while (i_locus<vec_full.size() && vec_full[i_locus]->bp < pos_bp[0]) {
		base = pos_cM[0];
		dist = (pos_bp[0] - vec_full[i_locus]->bp);
		vec_full[i_locus]->cm = base - mean_rate * dist;
		n_interpolated ++;
		i_locus ++;
	}

	//Set up middle positions using interpolation
	int closest_pos = 1;
	for (; i_locus < vec_full.size() ; ) {
		if (vec_full[i_locus]->cm == -1) {

			//Find suitable interpolation interval
			while (vec_full[i_locus]->bp > pos_bp[closest_pos] && closest_pos < pos_bp.size()) closest_pos++;

			//Interpolate
			if (closest_pos < pos_bp.size()) {
				assert(vec_full[i_locus]->bp < pos_bp[closest_pos]);
				assert(vec_full[i_locus]->bp > pos_bp[closest_pos-1]);
				base = pos_cM[closest_pos-1];
				rate = (pos_cM[closest_pos] - pos_cM[closest_pos-1]) / (pos_bp[closest_pos] - pos_bp[closest_pos-1]);
				dist = (vec_full[i_locus]->bp - pos_bp[closest_pos-1]);
				vec_full[i_locus]->cm = base + rate * dist;
				n_interpolated ++;
				i_locus ++;
			} else break;
		} else i_locus ++;
	}

	//Set up last positions to be mean rate
	while (i_locus < vec_full.size()) {
		base = pos_cM.back();
		dist = (vec_full[i_locus]->bp - pos_bp.back());
		vec_full[i_locus]->cm = base + mean_rate * dist;
		n_interpolated ++;
		i_locus ++;
	}
	return n_interpolated;
}

unsigned int variant_map::lengthBP() {
	return vec_full.back()->bp - vec_full[0]->bp + 1;
}

double variant_map::lengthcM() {
	return vec_full.back()->cm - vec_full[0]->cm;
}

void variant_map::setGeneticMap(gmap_reader & readerGM) {
	tac.clock();
	int n_set = setCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	int n_interpolated = interpolateCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	double baseline = vec_full[0]->cm;
	for (int l = 0 ; l < vec_full.size() ; l ++) vec_full[l]->cm -= baseline;
	vrb.bullet("cM interpolation [s=" + stb.str(n_set) + " / i=" + stb.str(n_interpolated) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("Region length [" + stb.str(vec_full.back()->bp-vec_full[0]->bp+1) + " bp / " + stb.str(vec_full.back()->cm-vec_full[0]->cm, 1) + " cM]");
}

void variant_map::setGeneticMap() {
	for (int l = 0 ; l < vec_full.size() ; l ++) vec_full[l]->cm = vec_full[l]->bp * 1.0 / 1e6;
	double baseline = vec_full[0]->cm;
	for (int l = 0 ; l < vec_full.size() ; l ++) vec_full[l]->cm -= baseline;
	vrb.bullet("Region length [" + stb.str(vec_full.back()->bp-vec_full[0]->bp+1) + " bp / " + stb.str(vec_full.back()->cm-vec_full[0]->cm, 1) + " cM (assuming 1cM per Mb)]");
}
