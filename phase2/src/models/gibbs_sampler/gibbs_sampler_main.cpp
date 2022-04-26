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

gibbs_sampler::gibbs_sampler(unsigned int _nsamples, unsigned int _niterations, unsigned int _nburnin) {
	nsamples = _nsamples;
	niterations = _niterations;
	nburnin = _nburnin;
	alleles = vector < bool > (2*nsamples, false);
	missing = vector < bool > (nsamples, false);
	phased = vector < bool > (nsamples, false);
	cstates = vector < vector < unsigned int > > (2*nsamples);
	cprobs = vector < vector < float > > (2*nsamples);
	pprobs = vector < float > (4*nsamples);
}

gibbs_sampler::~gibbs_sampler() {
	alleles.clear();
	missing.clear();
	unphased.clear();
	cstates.clear();
	cprobs.clear();
	pprobs.clear();
}

void gibbs_sampler::setHQ(){
	ee = 0.9999f;
	ed = 0.0001f;
}

void gibbs_sampler::setLQ(){
	ee = 0.9999f;
	ed = 0.0001f;
}

void gibbs_sampler::randomize_phase() {
	for(int h = 0 ; h < alleles.size() ; h +=2) {
		if (alleles[h+0] != alleles[h+1]) {
			if (rng.flipCoin()) {
				alleles[h+0] = !alleles[h+0];
				alleles[h+1] = !alleles[h+1];
			}
		}
	}
}

unsigned int fimax(vector < float > & vec) {
	float maxValue = vec[0];
	int maxIndex = 0;
	for (unsigned int i = 1; i < vec.size() ; i ++) if (vec[i] > maxValue) {
		maxValue = vec[i];
		maxIndex = i;
	}

	float epsilon = 1e-5;
	vector < int > max_indexes;
	for (unsigned int i = 0; i < vec.size() ; i ++) {
		if (abs(vec[i] - maxValue) < epsilon) {
			max_indexes.push_back(i);
		}
	}

	return max_indexes[rng.getInt(max_indexes.size())];
}

void gibbs_sampler::pushRare(genotype_set & G, unsigned int vr, unsigned int & n_yphased, unsigned int & n_nphased, float threshold) {
	vector < float > gprobs = vector < float > (4, 0.0f);
	for (int r = 0 ; r < G.GRvar_genotypes[vr].size() ; r ++) {
		if (G.GRvar_genotypes[vr][r].mis + G.GRvar_genotypes[vr][r].het) {
			if (!G.GRvar_genotypes[vr][r].pha) {
				copy(pprobs.begin() + 4 * G.GRvar_genotypes[vr][r].idx, pprobs.begin() + 4 * (G.GRvar_genotypes[vr][r].idx+1), gprobs.begin());
				int imax = fimax(gprobs);
				if (gprobs[imax] >= threshold) {
					switch (imax) {
					case 0: G.GRvar_genotypes[vr][r].al0 = 0; G.GRvar_genotypes[vr][r].al1 = 0; break;
					case 1: G.GRvar_genotypes[vr][r].al0 = 0; G.GRvar_genotypes[vr][r].al1 = 1; break;
					case 2: G.GRvar_genotypes[vr][r].al0 = 1; G.GRvar_genotypes[vr][r].al1 = 0; break;
					case 3: G.GRvar_genotypes[vr][r].al0 = 1; G.GRvar_genotypes[vr][r].al1 = 1; break;
					}
					G.GRvar_genotypes[vr][r].pha = 1;
					n_yphased ++;
				} else {
					n_nphased++;
				}
			} else n_yphased ++;
		}
	}
}

void gibbs_sampler::pushCommon(genotype_set & G, unsigned int vc, unsigned int & n_yphased, unsigned int & n_nphased, float threshold) {
	vector < float > gprobs = vector < float > (4, 0.0f);
	for (int i = 0 ; i < G.n_samples ; i ++) {
		if (missing[i] || (alleles[2*i+0] != alleles[2*i+1])) {
			if (!phased[i]) {
				copy(pprobs.begin() + 4 * i, pprobs.begin() + 4 * (i+1), gprobs.begin());
				int imax = fimax(gprobs);
				if (gprobs[imax] >= threshold) {
					switch (imax) {
					case 0: G.GCvar_alleles.set(vc, 2*i+0, false); G.GCvar_alleles.set(vc, 2*i+1, false); break;
					case 1: G.GCvar_alleles.set(vc, 2*i+0, false); G.GCvar_alleles.set(vc, 2*i+1, true); break;
					case 2: G.GCvar_alleles.set(vc, 2*i+0, true); G.GCvar_alleles.set(vc, 2*i+1, false); break;
					case 3: G.GCvar_alleles.set(vc, 2*i+0, true); G.GCvar_alleles.set(vc, 2*i+1, true); break;
					}
					G.GCvar_phased.set(vc, i, true);
					n_yphased ++;
				} else {
					n_nphased++;
				}
			} else n_yphased ++;
		}
	}
}

void gibbs_sampler::loadCommon(genotype_set & G, conditioning_set & C, state_set & P, unsigned int vc, float weight, bool use_compressed) {
	unsigned int vs = G.MAPC[vc];
	common = true;

	//Clean-up
	unphased.clear();
	for (int h = 0 ; h < C.n_haplotypes ; h ++) {
		cstates[h].clear();
		cprobs[h].clear();
	}
	fill (pprobs.begin(), pprobs.end() , 0.0f);

	//Genotype data
	for (int i = 0 ; i < G.n_samples ; i ++) {
		missing[i] = G.GCvar_missing.get(vc, i);
		phased[i] = G.GCvar_phased.get(vc, i);
		alleles[2*i+0] = G.GCvar_alleles.get(vc, 2*i+0);
		alleles[2*i+1] = G.GCvar_alleles.get(vc, 2*i+1);
		if (!phased[i]) {
			if (missing[i]) {
				alleles[2*i+0] = rng.flipCoin()?true:false;
				alleles[2*i+1] = rng.flipCoin()?true:false;
				unphased.push_back(i);
			} else if (alleles[2*i+0] != alleles[2*i+1]) {
				bool fc = rng.flipCoin();
				alleles[2*i+0] = fc?true:false;
				alleles[2*i+1] = fc?false:true;
				unphased.push_back(i);
			}
		}
	}

	//State probability data
	long int ci = P.Pmapping[vs];
	for (int u = 0 ; u < unphased.size() ; u ++) {
		for (int h = 0 ; h < 2 ; h ++) {
			//Iterate until hap is found

			if (use_compressed) {
				while (ci < P.Pstates1.size() && P.Pstates1[ci].id1 < (2*unphased[u]+h)) ci ++;
				assert(P.Pstates1[ci].id1 == (2*unphased[u]+h));
				for (int k = 0 ; ((ci+k) < P.Pstates1.size()) && (P.Pstates1[ci+k].id1 == (2*unphased[u]+h)) ; k++) {
					float pl = ((P.Pstates1[ci+k].lpb+1) * 1.0f) / 255;
					float pr = ((P.Pstates1[ci+k].rpb+1) * 1.0f) / 255;
					float kprob = pl * (1.0f - weight) + pr * weight;
					unsigned int kidx = C.indexes_pbwt_neighbour[2*unphased[u]+h][P.Pstates1[ci+k].kst];
					cstates[2*unphased[u]+h].push_back(kidx);
					cprobs[2*unphased[u]+h].push_back(kprob);
				}
			} else {
				while (ci < P.Pstates2.size() && P.Pstates2[ci].id1 < (2*unphased[u]+h)) ci ++;
				assert(P.Pstates2[ci].id1 == (2*unphased[u]+h));
				for (int k = 0 ; k < C.indexes_pbwt_neighbour[2*unphased[u]+h].size() ; k ++) {
					float pl = P.Pstates2[ci+k].lpb;
					float pr = P.Pstates2[ci+k].rpb;
					float kprob = pl * (1.0f - weight) + pr * weight;
					unsigned int kidx = C.indexes_pbwt_neighbour[2*unphased[u]+h][k];
					cstates[2*unphased[u]+h].push_back(kidx);
					cprobs[2*unphased[u]+h].push_back(kprob);
				}
			}
		}
	}
}

void gibbs_sampler::loadRare(genotype_set & G, conditioning_set & C, state_set & P, unsigned int vr, float weight, bool use_compressed) {
	unsigned int vs = G.MAPR[vr];
	common = false;

	//Clean-up
	unphased.clear();
	for (int h = 0 ; h < C.n_haplotypes ; h ++) {
		cstates[h].clear();
		cprobs[h].clear();
	}
	fill (pprobs.begin(), pprobs.end() , 0.0f);

	//Genotype data
	fill(alleles.begin(), alleles.end(), G.major_alleles[vr]);
	fill(missing.begin(), missing.end(), false);
	fill(phased.begin(), phased.end(), false);
	for (int r = 0 ; r < G.GRvar_genotypes[vr].size() ; r ++) {
		unsigned int ind = G.GRvar_genotypes[vr][r].idx;
		phased[ind] = G.GRvar_genotypes[vr][r].pha;

		if (phased[ind]) {
			alleles[2*ind+0] = G.GRvar_genotypes[vr][r].al0;
			alleles[2*ind+1] = G.GRvar_genotypes[vr][r].al1;
		} else {
			if (G.GRvar_genotypes[vr][r].mis) {
				missing[ind] = true;
				unphased.push_back(ind);
			} else if (G.GRvar_genotypes[vr][r].het) {
				bool fc = rng.flipCoin();
				alleles[2*ind+0] = fc?true:false;
				alleles[2*ind+1] = fc?false:true;
				unphased.push_back(ind);
			} else {
				alleles[2*ind+0] = !G.major_alleles[vr];
				alleles[2*ind+1] = !G.major_alleles[vr];
			}
		}
	}

	//State probability data
	long int ci = P.Pmapping[vs];
	for (int u = 0 ; u < unphased.size() ; u ++) {
		for (int h = 0 ; h < 2 ; h ++) {
			if (use_compressed) {
				while (ci < P.Pstates1.size() && P.Pstates1[ci].id1 < (2*unphased[u]+h)) ci ++;
				assert(P.Pstates1[ci].id1 == (2*unphased[u]+h));
				for (int k = 0 ; (ci+k) < P.Pstates1.size() && P.Pstates1[ci+k].id1 == (2*unphased[u]+h) ; k++) {
					float pl = ((P.Pstates1[ci+k].lpb+1) * 1.0f) / 255;
					float pr = ((P.Pstates1[ci+k].rpb+1) * 1.0f) / 255;
					float kprob = pl * (1.0f - weight) + pr * weight;
					unsigned int kidx = C.indexes_pbwt_neighbour[2*unphased[u]+h][P.Pstates1[ci+k].kst];
					cstates[2*unphased[u]+h].push_back(kidx);
					cprobs[2*unphased[u]+h].push_back(kprob);
				}
			} else {
				while (ci < P.Pstates2.size() && P.Pstates2[ci].id1 < (2*unphased[u]+h)) ci ++;
				assert(P.Pstates2[ci].id1 == (2*unphased[u]+h));
				for (int k = 0 ; k < C.indexes_pbwt_neighbour[2*unphased[u]+h].size() ; k ++) {
					float pl = P.Pstates2[ci+k].lpb;
					float pr = P.Pstates2[ci+k].rpb;
					float kprob = pl * (1.0f - weight) + pr * weight;
					unsigned int kidx = C.indexes_pbwt_neighbour[2*unphased[u]+h][k];
					cstates[2*unphased[u]+h].push_back(kidx);
					cprobs[2*unphased[u]+h].push_back(kprob);
				}
			}
		}
	}
}
