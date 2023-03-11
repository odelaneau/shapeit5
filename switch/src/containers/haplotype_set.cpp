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

#include <containers/haplotype_set.h>

using namespace std;

haplotype_set::haplotype_set() {
	clear();
}

haplotype_set::~haplotype_set() {
	clear();
}

void haplotype_set::clear() {
	n_variants = 0;
	Htrue.clear();
	Hesti.clear();
	IDXesti.clear();
	Missing.clear();
	Phased.clear();
	MAC.clear();
	Mothers.clear();
	Fathers.clear();
}

void haplotype_set::push(string & sample_id) {

	Htrue.push_back(vector < bool > ());
	Htrue.push_back(vector < bool > ());
	Hesti.push_back(vector < bool > ());
	Hesti.push_back(vector < bool > ());
	Hprob.push_back(vector < bool > ());
	Missing.push_back(vector < bool > ());
	Phased.push_back(vector < bool > ());
	Estimated.push_back(vector < bool > ());

	mapSamples.insert(pair < string, int > (sample_id, vecSamples.size()));
	vecSamples.push_back(sample_id);
	Mothers.push_back(-1);
	Fathers.push_back(-1);
}


void haplotype_set::readPedigrees(string fped, bool dupid) {
	string buffer;
	vector < string > str;
	vrb.title("Read pedigrees in [" + fped + "]");
	input_file fd (fped);
	if (fd.fail()) vrb.error("Cannot open file!");
	int n_unr = 0, n_duo = 0, n_tri = 0;
	while (getline(fd, buffer)) {
		stb.split(buffer, str);
		if (dupid) for (int i = 0 ; i < 3 ; i++) str[i] = str[i] + "_" + str[i];
		map < string, int > :: iterator itC = mapSamples.find(str[0]);
		map < string, int > :: iterator itF = mapSamples.find(str[1]);
		map < string, int > :: iterator itM = mapSamples.find(str[2]);
		if (itC != mapSamples.end()) {
			if (itF != mapSamples.end()) Fathers[itC->second] = itF->second;
			if (itM != mapSamples.end()) Mothers[itC->second] = itM->second;
			int type = (Fathers[itC->second] >= 0) + (Mothers[itC->second] >= 0);
			switch (type) {
				case 2: n_tri ++; break;
				case 1: n_duo ++; break;
				default: n_unr ++; break;
			}
		}
	}
	fd.close();
	vrb.bullet("#trios = " + stb.str(n_tri));
	vrb.bullet("#duos = " + stb.str(n_duo));
	vrb.bullet("#unrelateds = " + stb.str(n_unr));
}

void haplotype_set::assumePhased() {
	vrb.title("Assuming all hets in validation are correctly phased");
	Phased = vector < vector < bool > > (vecSamples.size(), vector < bool > (n_variants, false));
	for (int i = 0 ; i < vecSamples.size() ; i ++) for (int l = 0 ; l < n_variants ; l ++) Phased[i][l] = !Missing[i][l] && (Htrue[2*i+0][l]!=Htrue[2*i+1][l]);
}

unsigned int haplotype_set::distance(unsigned int h0, unsigned int h1) {
	unsigned int dist = 0;
	for (int l = 0 ; l < n_variants ; l++) dist += (Phased[h1/2][l] && (Htrue[h0][l] != Hesti[h1][l]));
	return dist;
}

