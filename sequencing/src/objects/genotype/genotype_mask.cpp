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
#include <objects/genotype/genotype_header.h>

void genotype::mask() {
	// Check if there is PS information
	bool toBeProcessed = false;
	ProbabilityMask.clear();
	for (int p = 0 ; p < PhaseSets.size() ; p ++ ) if (PhaseSets[p].ps > 0) {
		toBeProcessed=true;
	}
	if (toBeProcessed) {
		// Allocate ProbabilityMask
		ProbabilityMask = vector < bool > (n_transitions, true);
		// Iterates over segments
		unsigned char prev_dipcodes [64];
		unsigned int prev_dipcount = 1, curr_dipcount = 0, n_curr_trans = 0, n_curr_amb = 0;
		for (unsigned int s = 0, a = 0, v = 0, t = 0 ; s < n_segments ; s ++) {

			if (s == 0) {
				// Get numbers of a, v, and t
				curr_dipcount = countDiplotypes(Diplotypes[s]);
				makeDiplotypes(Diplotypes[s]);
				n_curr_trans = curr_dipcount * prev_dipcount;
				for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++) n_curr_amb += (VAR_GET_AMB(MOD2(v+vrel), Variants[DIV2(v+vrel)]));

				// Build haplotype blocks from PS informations
				map < int, vector < char > > HapBlocks;
				for (unsigned int vrel = 0, arel = 0 ; vrel < Lengths[s] ; vrel ++) {
					if (VAR_GET_AMB(MOD2(v+vrel), Variants[DIV2(v+vrel)])) {
						unsigned int ps = PhaseSets[a+arel].ps;
						if (ps) {
							bool allele0 = PhaseSets[a+arel].a0;
							bool allele1 = PhaseSets[a+arel].a1;
							map < int, vector < char > >::iterator it = HapBlocks.find(ps);
							if (it != HapBlocks.end()) {
								it->second[2*arel+0] = allele0;
								it->second[2*arel+1] = allele1;
							} else {
								vector < char > newHapBlock = vector < char > (2*n_curr_amb, -1);
								newHapBlock[2*arel+0] = allele0;
								newHapBlock[2*arel+1] = allele1;
								HapBlocks.insert(pair < int, vector < char > > (ps, newHapBlock));
							}
						}
						arel ++;
					}
				}
/*
				cout << "BLOC:" << endl;
				for (map < int , vector < char > > :: iterator it = HapBlocks.begin() ; it != HapBlocks.end() ; ++it) {
					for (int l = 0 ; l < it->second.size() ; l += 2) cout << "\t" << (int)it->second[l];
					cout << endl;
					for (int l = 0 ; l < it->second.size() ; l += 2) cout << "\t" << (int)it->second[l+1];
					cout << endl;
				}
*/
				// Now, iterates over transitions and check consistencies
				vector < char > haplotype = vector < char > (n_curr_amb, 0);
				for (int trel = 0 ; trel < n_curr_trans ; trel ++) {
					unsigned int curr_dip = curr_dipcodes[trel];
					unsigned int curr_h0 = DIP_HAP0(curr_dip);
					for (unsigned int vrel = 0, arel = 0 ; vrel < Lengths[s] ; vrel ++) {
						if (VAR_GET_AMB(MOD2(v+vrel), Variants[DIV2(v+vrel)])) {
							haplotype[arel] = HAP_GET(Ambiguous[a+arel], curr_h0);
							arel ++;
						}
					}
					for (map < int , vector < char > > :: iterator it = HapBlocks.begin() ; it != HapBlocks.end() ; ++it) {
						assert(it->second.size() == 2 * haplotype.size());
						int first_allele = -1;
						for (unsigned int l = 0 ; l < n_curr_amb ; l ++) {
							if (it->second[2*l] >= 0) {
								if (first_allele < 0) first_allele = (haplotype[l] != it->second[2*l]);
								else if (haplotype[l] != it->second[2*l+first_allele]) ProbabilityMask[t+trel] = false;
							}
						}
					}
				}

/*
				double prop = 0.0;
				for (int trel = 0 ; trel < n_curr_trans ; trel ++) prop += ProbabilityMask[t+trel];
				cout << name << "\t" << s << "\t" << HapBlocks.size() << "\t" << prop << "\t" << prop * 1.0 / n_curr_trans << endl;
*/
				// Update t cursor
				t += n_curr_trans;
				std::copy(curr_dipcodes, curr_dipcodes + curr_dipcount, prev_dipcodes);
				prev_dipcount = curr_dipcount;

			} else {
				// Get numbers of a, v, and t
				n_curr_amb = 0;
				curr_dipcount = countDiplotypes(Diplotypes[s]);
				makeDiplotypes(Diplotypes[s]);
				n_curr_trans = curr_dipcount * prev_dipcount;
				for (unsigned int vrel = 0 ; vrel < Lengths[s-1] + Lengths[s] ; vrel ++) n_curr_amb += (VAR_GET_AMB(MOD2(v+vrel), Variants[DIV2(v+vrel)]));

				// Build haplotype blocks from PS informations
				map < int, vector < char > > HapBlocks;
				for (unsigned int vrel = 0, arel = 0 ; vrel < Lengths[s-1] + Lengths[s] ; vrel ++) {
					if (VAR_GET_AMB(MOD2(v+vrel), Variants[DIV2(v+vrel)])) {
						unsigned int ps = PhaseSets[a+arel].ps;
						if (ps) {
							bool allele0 = PhaseSets[a+arel].a0;
							bool allele1 = PhaseSets[a+arel].a1;
							map < int, vector < char > >::iterator it = HapBlocks.find(ps);
							if (it != HapBlocks.end()) {
								it->second[2*arel+0] = allele0;
								it->second[2*arel+1] = allele1;
							} else {
								vector < char > newHapBlock = vector < char > (2*n_curr_amb, -1);
								newHapBlock[2*arel+0] = allele0;
								newHapBlock[2*arel+1] = allele1;
								HapBlocks.insert(pair < int, vector < char > > (ps, newHapBlock));
							}
						}
						arel ++;
					}
				}
/*
				cout << "BLOC:" << endl;
				for (map < int , vector < char > > :: iterator it = HapBlocks.begin() ; it != HapBlocks.end() ; ++it) {
					for (int l = 0 ; l < it->second.size() ; l += 2) cout << "\t" << (int)it->second[l];
					cout << endl;
					for (int l = 0 ; l < it->second.size() ; l += 2) cout << "\t" << (int)it->second[l+1];
					cout << endl;
				}
*/
				// Now, iterates over transitions and check consistencies
				vector < char > haplotype = vector < char > (n_curr_amb, 0);
				for (int trel = 0 ; trel < n_curr_trans ; trel ++) {
					unsigned int prev_dip = prev_dipcodes[trel/curr_dipcount];
					unsigned int next_dip = curr_dipcodes[trel%curr_dipcount];
					unsigned int prev_h0 = DIP_HAP0(prev_dip);
					unsigned int next_h0 = DIP_HAP0(next_dip);
					for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++) {
						if (VAR_GET_AMB(MOD2(v+vrel), Variants[DIV2(v+vrel)])) {
							haplotype[arel] = HAP_GET(Ambiguous[a+arel], (vrel<Lengths[s-1])?prev_h0:next_h0);
							arel ++;
						}
					}
					for (map < int , vector < char > > :: iterator it = HapBlocks.begin() ; it != HapBlocks.end() ; ++it) {
						assert(it->second.size() == 2 * haplotype.size());
						int first_allele = -1;
						for (unsigned int l = 0 ; l < n_curr_amb ; l ++) {
							if (it->second[2*l] >= 0) {
								if (first_allele < 0) first_allele = (haplotype[l] != it->second[2*l]);
								else if (haplotype[l] != it->second[2*l+first_allele]) ProbabilityMask[t+trel] = false;
							}
						}
					}
				}
/*
				double prop = 0.0;
				for (int trel = 0 ; trel < n_curr_trans ; trel ++) prop += ProbabilityMask[t+trel];
				cout << name << "\t" << s << "\t" << HapBlocks.size() << "\t" << prop << "\t" << prop * 1.0 / n_curr_trans << endl;
*/
				// Update a, v and t cursors
				for (unsigned int vrel = 0 ; vrel < Lengths[s-1] ; vrel ++) a += VAR_GET_AMB(MOD2(v+vrel), Variants[DIV2(v+vrel)]);
				v += Lengths[s-1];
				t += n_curr_trans;
				std::copy(curr_dipcodes, curr_dipcodes + curr_dipcount, prev_dipcodes);
				prev_dipcount = curr_dipcount;
			}
		}
	}
}
