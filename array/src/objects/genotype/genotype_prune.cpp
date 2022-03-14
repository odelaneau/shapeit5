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

class Transition {
public:
	double prob;
	unsigned int idx;
	Transition() { prob = 0.0; idx = 0;}
	~Transition() {}
	bool operator < (const Transition & t) const { return prob > t.prob; }
};

class TransStatistics {
public:
	double entropy;
	unsigned int idx;
	bool merged;
	TransStatistics() {entropy = 1000; idx = -1; merged = false; }
	~TransStatistics() {};
	bool operator < (const TransStatistics & s) const { return entropy < s.entropy; }
};

void genotype::mapMerges(vector < double > & currProbs, double thresholdProbMass, vector < bool > & flagMerges) {
	vector < TransStatistics > vecTransStatistics = vector < TransStatistics > (n_segments - 1);
	vector < Transition > vecTransitions = vector < Transition > (4096);

	//Step0: initialize cursors
	unsigned int prev_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned int curr_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned char prev_dipcodes [64];
	makeDiplotypes(Diplotypes[0]);
	std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
	unsigned int toffset = prev_dipcount;
	unsigned int n_curr_transitions = 0;
	unsigned int aoffset = 0, voffset = 0;

	for (int s = 1 ; s < n_segments ; s++) {
		//Step1: update cursors (1)
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		n_curr_transitions = prev_dipcount * curr_dipcount;
		makeDiplotypes(Diplotypes[s]);

		//Step2: intialize transition statistics
		vecTransStatistics[s-1].idx = s;
		vecTransStatistics[s-1].merged = false;
		vecTransStatistics[s-1].entropy = 4096;

		//Step3: check number of variants in merged segment
		unsigned int segment_length = Lengths[s-1] + Lengths[s];
		if (segment_length < std::numeric_limits< unsigned short >::max()) {
			unsigned int n_ambiguous_merged = 0;
			for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++)
				if (VAR_GET_AMB(MOD2(voffset+vrel), Variants[DIV2(voffset+vrel)]))
					n_ambiguous_merged++;
			//Step4: check number of ambiguous variants in merged segment
			if (n_ambiguous_merged < MAX_AMB) {
				//Step5: sort transitions by decreasing order
				for (int t = 0 ; t < n_curr_transitions ; t ++) {
					vecTransitions[t].prob = currProbs[toffset + t];
					vecTransitions[t].idx = t;
				}
				sort(vecTransitions.begin(), vecTransitions.begin() + n_curr_transitions);
				//Step6: compute transition entropy
				vecTransStatistics[s-1].entropy = 0.0;
				for (int t = 0 ; t < n_curr_transitions ; t ++) {
					double cProb = vecTransitions[t].prob;
					double lProb = -1.0 * ((cProb==0.0f)?0:log10(cProb));
					vecTransStatistics[s-1].entropy += cProb * lProb;
				}
				//Step7: check that 8 haplotypes capture lots of the cumulative probability mass
				double cumSumProbs = 0.0;
				vector < int > Mhaps = vector < int >(HAP_NUMBER * HAP_NUMBER, -1);
				for (int t = 0, n_haps = 0 ; t < n_curr_transitions ; t ++) {
					cumSumProbs += vecTransitions[t].prob;
					unsigned int prev_dip = prev_dipcodes[vecTransitions[t].idx/curr_dipcount];
					unsigned int next_dip = curr_dipcodes[vecTransitions[t].idx%curr_dipcount];
					unsigned int prev_h0 = DIP_HAP0(prev_dip);
					unsigned int prev_h1 = DIP_HAP1(prev_dip);
					unsigned int next_h0 = DIP_HAP0(next_dip);
					unsigned int next_h1 = DIP_HAP1(next_dip);
					unsigned int merged_h0 = prev_h0 * HAP_NUMBER + next_h0;
					unsigned int merged_h1 = prev_h1 * HAP_NUMBER + next_h1;
					bool new_h0 = (Mhaps[merged_h0] < 0);
					bool new_h1 = ((Mhaps[merged_h1] < 0) && (merged_h0 != merged_h1));
					//if ((n_haps + new_h0 + new_h1) <= HAP_NUMBER) {
						if (new_h0) Mhaps[merged_h0] = n_haps++;
						if (new_h1) Mhaps[merged_h1] = n_haps++;
					//}
					if (n_haps == HAP_NUMBER && cumSumProbs > thresholdProbMass) vecTransStatistics[s-1].merged = true;
				}
			}
		}

		//Step8: update cursors (2)
		for (unsigned int vrel = 0 ; vrel < Lengths[s-1] ; vrel ++) if (VAR_GET_AMB(MOD2(voffset+vrel), Variants[DIV2(voffset+vrel)])) aoffset++;
		voffset += Lengths[s-1];
		std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
		prev_dipcount = curr_dipcount;
		toffset += n_curr_transitions;
	}
	//Step9: map acceptable merges
	sort(vecTransStatistics.begin(), vecTransStatistics.end());
	flagMerges = vector < bool > (n_segments+1, false);
	for (unsigned int s = 0 ; s < vecTransStatistics.size() ; s ++) {
		bool no_adjacent_merges = !flagMerges[vecTransStatistics[s].idx-1] && !flagMerges[vecTransStatistics[s].idx+1];
		bool can_be_merged = vecTransStatistics[s].merged;
		flagMerges[vecTransStatistics[s].idx] = (no_adjacent_merges && can_be_merged);
	}
}

void genotype::performMerges(vector < double > & currProbs, vector < bool > & flagMerges) {
	vector < Transition > vecTransitions = vector < Transition > (4096);

	//Step0: initialize duplicates
	vector < unsigned char > Ambiguous2 = vector < unsigned char > (Ambiguous.size(), 0);
	vector < unsigned long > Diplotypes2;
	vector < unsigned short > Lengths2;
	unsigned int n_segments2 = n_segments;
	for (int s = 0 ; s < flagMerges.size() ; s++) n_segments2 -= flagMerges[s];
	Diplotypes2.reserve(n_segments2);
	Lengths2.reserve(n_segments2);

	//Step1: initialize cursors
	unsigned int prev_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned int curr_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned char prev_dipcodes [64];
	makeDiplotypes(Diplotypes[0]);
	std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
	unsigned int toffset = prev_dipcount;
	unsigned int n_curr_transitions = 0;
	unsigned int aoffset = 0, voffset = 0;

	for (int s = 1 ; s < flagMerges.size() -1 ; s ++) {
		//Step1: update cursors (1)
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		n_curr_transitions = prev_dipcount * curr_dipcount;
		makeDiplotypes(Diplotypes[s]);

		//case1: merge to be done
		if (flagMerges[s]) {
			//cout << name << " M " << aoffset << endl;
			Lengths2.push_back(Lengths[s-1]+Lengths[s]);
			Diplotypes2.push_back(0x0000000000000000UL);
			for (int t = 0 ; t < n_curr_transitions ; t ++) { vecTransitions[t].prob = currProbs[toffset + t]; vecTransitions[t].idx = t; }
			sort(vecTransitions.begin(), vecTransitions.begin() + n_curr_transitions);
			int n_haps = 0;
			vector < int > Mhaps = vector < int >(HAP_NUMBER * HAP_NUMBER, -1);
			for (int t = 0 ; t < n_curr_transitions ; t ++) {
				unsigned int prev_dip = prev_dipcodes[vecTransitions[t].idx/curr_dipcount];
				unsigned int next_dip = curr_dipcodes[vecTransitions[t].idx%curr_dipcount];
				unsigned int prev_h0 = DIP_HAP0(prev_dip);
				unsigned int prev_h1 = DIP_HAP1(prev_dip);
				unsigned int next_h0 = DIP_HAP0(next_dip);
				unsigned int next_h1 = DIP_HAP1(next_dip);
				unsigned int merged_h0 = prev_h0 * HAP_NUMBER + next_h0;
				unsigned int merged_h1 = prev_h1 * HAP_NUMBER + next_h1;
				bool new_h0 = (Mhaps[merged_h0] < 0);
				bool new_h1 = ((Mhaps[merged_h1] < 0) && (merged_h0 != merged_h1));
				if ((n_haps + new_h0 + new_h1) <= HAP_NUMBER) {
					if (new_h0) {
						Mhaps[merged_h0] = n_haps;
						for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++) {
							if (VAR_GET_AMB(MOD2(voffset+vrel), Variants[DIV2(voffset+vrel)])) {
								if (HAP_GET(Ambiguous[aoffset+arel], (vrel<Lengths[s-1])?prev_h0:next_h0)) HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h0]);
								arel ++;
							}
						}
						n_haps ++;
					}
					if (new_h1) {
						Mhaps[merged_h1] = n_haps;
						for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++) {
							if (VAR_GET_AMB(MOD2(voffset+vrel), Variants[DIV2(voffset+vrel)])) {
								if (HAP_GET(Ambiguous[aoffset+arel], (vrel<Lengths[s-1])?prev_h1:next_h1)) HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h1]);
								arel ++;
							}
						}
						n_haps ++;
					}
					DIP_SET(Diplotypes2.back(), Mhaps[merged_h0] * HAP_NUMBER + Mhaps[merged_h1]);
				}
			}
			assert(n_haps == HAP_NUMBER);
		//Case2: no merge to be done, push last segment
		} else if (!flagMerges[s-1]) {
			//cout << name << " C " << aoffset << endl;
			for (unsigned int vrel = 0, arel = 0 ; vrel < Lengths[s-1] ; vrel ++) {
				if (VAR_GET_AMB(MOD2(voffset+vrel), Variants[DIV2(voffset+vrel)])) {
					Ambiguous2[aoffset+arel] = Ambiguous[aoffset+arel];
					arel ++;
				}
			}
			Lengths2.push_back(Lengths[s-1]);
			Diplotypes2.push_back(Diplotypes[s-1]);
		}

		//Update cursors
		for (unsigned int vrel = 0 ; vrel < Lengths[s-1] ; vrel ++) if (VAR_GET_AMB(MOD2(voffset+vrel), Variants[DIV2(voffset+vrel)])) aoffset++;
		voffset += Lengths[s-1];
		std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
		prev_dipcount = curr_dipcount;
		toffset += n_curr_transitions;
	}
	if (!flagMerges[flagMerges.size()-2]) {
		for (unsigned int vrel = 0, arel = 0 ; vrel < Lengths.back() ; vrel ++) {
			if (VAR_GET_AMB(MOD2(voffset+vrel), Variants[DIV2(voffset+vrel)])) {
				Ambiguous2[aoffset+arel] = Ambiguous[aoffset+arel];
				arel ++;
			}
		}
		Lengths2.push_back(Lengths.back());
		Diplotypes2.push_back(Diplotypes.back());
	}

	//free();
	Ambiguous = Ambiguous2;
	Diplotypes = Diplotypes2;
	Lengths = Lengths2;
	n_segments = n_segments2;
	n_transitions = countTransitions();
}
