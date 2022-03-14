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

void genotype::sample(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {
	if (rng.getDouble() < 0.5f) sampleForward(CurrentTransProbabilities, CurrentMissingProbabilities);
	else sampleBackward(CurrentTransProbabilities, CurrentMissingProbabilities);
}

void genotype::sampleForward(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {
	double sumProbs = 0.0;
	unsigned int prev_sampled = 0;
	unsigned int curr_dipcount = 0, prev_dipcount = 1;
	vector < double > currProbs = vector < double > (64, 0.0);
	vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);
	for (unsigned int s = 0, toffset = 0 ; s < n_segments ; s ++) {
		sumProbs = 0.0;
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		for (unsigned int tabs = toffset + prev_sampled*curr_dipcount, trel = 0 ; trel < curr_dipcount ; ++trel, ++tabs)
			sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
		prev_sampled = rng.sample(currProbs, sumProbs);
		makeDiplotypes(Diplotypes[s]);
		DipSampled[s] = curr_dipcodes[prev_sampled];
		toffset += prev_dipcount * curr_dipcount;
		prev_dipcount = curr_dipcount;
	}
	make(DipSampled, CurrentMissingProbabilities);
}

void genotype::sampleBackward(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {

	double sumProbs = 0.0;
	int next_sampled = -1;
	unsigned int curr_dipcount = 0, next_dipcount = countDiplotypes(Diplotypes[n_segments - 1]);
	vector < double > currProbs = vector < double > (64 * 64, 0.0);
	vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);

	for (int s = n_segments - 2, toffset = n_transitions ; s >= 0 ; s --) {
		sumProbs = 0.0;
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		toffset -= next_dipcount * curr_dipcount;

		if (next_sampled >= 0) {
			currProbs.resize(64);
			for (unsigned int tabs = toffset+next_sampled, trel = 0 ; trel < curr_dipcount ; ++trel, tabs += next_dipcount)
				sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
			next_sampled = rng.sample(currProbs, sumProbs);
			makeDiplotypes(Diplotypes[s]);
			DipSampled[s] = curr_dipcodes[next_sampled];
		} else {
			for (unsigned int tabs = toffset, trel = 0 ; tabs < n_transitions ; ++trel, ++tabs)
				sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
			next_sampled = rng.sample(currProbs, sumProbs);
			makeDiplotypes(Diplotypes[s+1]);
			DipSampled[s+1] = curr_dipcodes[next_sampled % next_dipcount];
			makeDiplotypes(Diplotypes[s]);
			next_sampled = next_sampled / next_dipcount;
			DipSampled[s] = curr_dipcodes[next_sampled];
		}
		next_dipcount = curr_dipcount;
	}
	make(DipSampled, CurrentMissingProbabilities);
}

void genotype::solve() {
	unsigned int curr_dipcount = 0, prev_dipcount = 1;
	vector < vector < double > > maxProbs = vector < vector < double > > (n_segments, vector < double > ());
	vector < vector < int > > maxIndexes = vector < vector < int > > (n_segments, vector < int > ());

	for (int s = 0, toffset = 0, trel = 0 ; s < n_segments ; s ++) {
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		maxProbs[s] = vector < double > (curr_dipcount, 0.0);
		maxIndexes[s] = vector < int > (curr_dipcount, 0);
		for (int t = 0 ; t < prev_dipcount * curr_dipcount ; t++) {
			int prev_dip = t/curr_dipcount;
			int next_dip = t%curr_dipcount;
			//double currProb = (s?maxProbs[s-1][prev_dip]:1.0) * StoredProbs[t+toffset];
			double currProb = (s?maxProbs[s-1][prev_dip]:1.0) * (ProbMask[t+toffset]?ProbStored[trel++]:1e-6);
			if (currProb > maxProbs[s][next_dip]) {
				maxProbs[s][next_dip] = currProb;
				maxIndexes[s][next_dip] = prev_dip;
			}
		}
		double sumProb = 0.0;
		for (int d = 0 ; d < curr_dipcount ; d ++) sumProb += maxProbs[s][d];
		for (int d = 0 ; d < curr_dipcount ; d ++) maxProbs[s][d] /= sumProb;
		toffset += prev_dipcount * curr_dipcount;
		prev_dipcount = curr_dipcount;
	}

	vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);
	unsigned int bestDip = alg.imax(maxProbs.back());
	makeDiplotypes(Diplotypes.back());
	DipSampled.back() = curr_dipcodes[bestDip];
	for (int s = DipSampled.size() - 2 ; s >= 0 ; s --) {
		bestDip = maxIndexes[s+1][bestDip];
		makeDiplotypes(Diplotypes[s]);
		DipSampled[s] = curr_dipcodes[bestDip];
	}
	make(DipSampled);
}

void genotype::store(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {
	if (ProbMask.size() == 0) {
		n_stored_transitionProbs = 0;
		ProbMask = vector < bool > (n_transitions, false);
		for (unsigned int t = 0 ; t < n_transitions ; t ++) if (CurrentTransProbabilities[t] >= 1e-6) {
			n_stored_transitionProbs ++;
			ProbMask[t] = true;
		}
		ProbStored = vector  < float > (n_stored_transitionProbs, 0.0);
		ProbMissing = vector < float > (n_missing * HAP_NUMBER, 0.0);
	}
	for (unsigned int t = 0, trel = 0 ; t < n_transitions ; t ++) {
		if (ProbMask[t]) ProbStored[trel++] += CurrentTransProbabilities[t];
		//if (!ProbMask[t] && CurrentTransProbabilities[t] > 0.01) cout << "Should be stored " << CurrentTransProbabilities[t] << endl;
	}
	for (unsigned int m = 0 ; m < (n_missing * HAP_NUMBER) ; m ++) ProbMissing[m] += CurrentMissingProbabilities[m];
	n_storage_events ++;
}
