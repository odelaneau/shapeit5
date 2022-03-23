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

#include <models/gibbs_sampler/gibbs_sampler_header.h>

void gibbs_sampler::iterate(int & error, int & total) {
	float sum;
	vector < float > hprob = vector < float > (4, 0.0f);
	vector < float > gprob = vector < float > (4, 0.0f);

	for (int iter = 0 ; iter < niterations ; iter ++) {

		//cout << " ========================================== " << endl;

		for (int ui = 0 ; ui < unphased.size() ; ui ++) {

			//Init
			int ind = unphased[ui];
			fill(hprob.begin(), hprob.end(), 0.0f);
			fill(gprob.begin(), gprob.end(), 1.0f);
			assert(missing [ind] || (alleles[2*ind+0] != alleles[2*ind+1]));

			//
			for (int h = 0 ; h < 2 ; h++) {
				assert(cstates[2*ind+h].size());
				for (int k = 0 ; k < cstates[2*ind+h].size() ; k ++) {
					bool a = alleles[cstates[2*ind+h][k]];
					float p = cprobs[2*ind+h][k];
					hprob[2*h+a] += p;
				}
				sum = hprob[2*h+0] + hprob[2*h+1];
				hprob[2*h+0] /= sum;
				hprob[2*h+1] /= sum;
			}

			for (int h = 0 ; h < 4 ; h++) if (hprob[h] < 1e-7) hprob[h] = 1e-7;

			//
			if (!missing [ind]) {
				gprob[0] = 0.0f;
				gprob[3] = 0.0f;
			}
			gprob[0] *= (hprob[0] * hprob[2]);
			gprob[1] *= (hprob[0] * hprob[3]);
			gprob[2] *= (hprob[1] * hprob[2]);
			gprob[3] *= (hprob[1] * hprob[3]);
			sum = accumulate(gprob.begin(), gprob.end(), 0.0f);
			//for (int g = 0 ; g < 4 ; g++) gprob[g] /= sum;

			//cout << iter << " " << ui << " " << stb.str(sum, 3) << " / " << stb.str(gprob, 3) << " / " << stb.str(hprob, 3) << endl;
			assert(sum > 0);

			//
			int sampleg = rng.sample(gprob, sum);
			switch (sampleg) {
			case 0:	alleles[2*ind+0] = 0; alleles[2*ind+1] = 0; break;
			case 1:	alleles[2*ind+0] = 0; alleles[2*ind+1] = 1; break;
			case 2:	alleles[2*ind+0] = 1; alleles[2*ind+1] = 0; break;
			case 3:	alleles[2*ind+0] = 1; alleles[2*ind+1] = 1; break;
			}

			//
			if (alleles[2*ind+0] == alleles[2*ind+1]) {
				cout << iter << " " << ui << " " << sampleg << " " << stb.str(sum, 3) << " / " << stb.str(gprob, 3) << " / " << stb.str(hprob, 3) << endl;
				exit(1);
			}

			//
			if (iter >= nburnin) {
				pprobs[4*ind+0] += gprob[0];
				pprobs[4*ind+1] += gprob[1];
				pprobs[4*ind+2] += gprob[2];
				pprobs[4*ind+3] += gprob[3];
			}
		}
	}
	for (int ui = 0 ; ui < unphased.size() ; ui ++) {
		int ind = unphased[ui];
		unsigned int bestg = distance(pprobs.begin() + 4*ind, max_element(pprobs.begin() + 4*ind, pprobs.begin() + 4*(ind+1))) % 4;
		switch (bestg) {
		case 0:	alleles[2*ind+0] = 0; alleles[2*ind+1] = 0; break;
		case 1:	alleles[2*ind+0] = 0; alleles[2*ind+1] = 1; break;
		case 2:	alleles[2*ind+0] = 1; alleles[2*ind+1] = 0; break;
		case 3:	alleles[2*ind+0] = 1; alleles[2*ind+1] = 1; break;
		}
		if (unphased.size() > 1) {
			error +=  (alleles[2*ind+0] != truth[2*ind+0]);
			assert(alleles[2*ind+0] != alleles[2*ind+1]);
			assert(truth[2*ind+0] != truth[2*ind+1]);
			total++;
		}
	}
}
