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

gibbs_sampler::gibbs_sampler() {
}

gibbs_sampler::~gibbs_sampler() {
}

void gibbs_sampler::allocate(unsigned int _nsamples, unsigned int _niterations, unsigned int _nburnin) {
	nsamples = _nsamples;
	niterations = _niterations;
	nburnin = _nburnin;
	alleles0 = vector < bool > (nsamples);
	alleles1 = vector < bool > (nsamples);
	missing = vector < bool > (nsamples);
	state_indexes = vector < vector < unsigned int > > (2 * nsamples);
	state_lprobs = vector < vector < float > > (2 * nsamples);
	state_rprobs = vector < vector < float > > (2 * nsamples);
	phasing_probs = vector < float > (4 * nsamples);
}

bool gibbs_sampler::loadCommonUnphasedGenotypes(unsigned int vc, genotype_set & GS) {
	unphased_indexes.clear();
	for (int i = 0 ; i < GS.n_samples ; i++) {
		missing[i] = GS.GCmissing[vc][i];
		alleles0[i] = GS.GCalleles[vc][2*i+0];
		alleles1[i] = GS.GCalleles[vc][2*i+1];
		if (missing[i]) {
			unphased_indexes.push_back(i);
			alleles0[i] = rng.flipCoin()?false:true;
			alleles1[i] = rng.flipCoin()?false:true;
		} else if (alleles0[i] != alleles1[i]) {
			unphased_indexes.push_back(i);
			bool rcoin = rng.flipCoin();
			alleles0[i] = rcoin?true:false;
			alleles1[i] = rcoin?false:true;
		}
	}
	return (unphased_indexes.size() > 0);
}

bool gibbs_sampler::loadRareUnphasedGenotypes(unsigned int vr, genotype_set & GS, bool major_allele) {
	unphased_indexes.clear();
	std::fill(alleles0.begin(), alleles0.end(), major_allele);
	std::fill(alleles1.begin(), alleles1.end(), major_allele);
	std::fill(missing.begin(), missing.end(), false);
	for (int i = 0 ; i < GS.GRindexes[vr].size() ; i++) {
		if (GS.GRmissing[vr][i]) {
			unphased_indexes.push_back(GS.GRindexes[vr][i]);
			missing[GS.GRindexes[vr][i]] = true;
			alleles0[GS.GRindexes[vr][i]] = major_allele;
			alleles1[GS.GRindexes[vr][i]] = major_allele;
		} else if (GS.GRhets[vr][i]) {
			unphased_indexes.push_back(GS.GRindexes[vr][i]);
			bool rcoin = rng.flipCoin();
			alleles0[GS.GRindexes[vr][i]] = rcoin?true:false;
			alleles1[GS.GRindexes[vr][i]] = rcoin?false:true;
		} else {
			alleles0[GS.GRindexes[vr][i]] = !major_allele;
			alleles1[GS.GRindexes[vr][i]] = !major_allele;
		}
	}
	return (unphased_indexes.size() > 0);
}

void gibbs_sampler::pushCommonPhasedGenotypes(unsigned int vc, genotype_set & GS) {
	for (int ui = 0 ; ui < unphased_indexes.size() ; ui++) {
		GS.GCalleles[vc][2*unphased_indexes[ui]+0] = alleles0[unphased_indexes[ui]];
		GS.GCalleles[vc][2*unphased_indexes[ui]+1] = alleles1[unphased_indexes[ui]];
	}
}

void gibbs_sampler::pushRarePhasedGenotypes(unsigned int vr, genotype_set & GS) {
	for (int i = 0 ; i < GS.GRindexes[vr].size() ; i++) {
		if (GS.GRmissing[vr][i] || GS.GRhets[vr][i]) {
			GS.GRalleles[vr][2*i+0] = alleles0[GS.GRindexes[vr][i]];
			GS.GRalleles[vr][2*i+1] = alleles1[GS.GRindexes[vr][i]];
		}
	}
}

void gibbs_sampler::loadStateSpace(unsigned int hap, vector < unsigned int > & states, vector < float > & lprod, vector < float > & rprod, float threshold) {
	state_indexes[hap].clear();
	state_lprobs[hap].clear();
	state_rprobs[hap].clear();
	for (int k = 0 ; k < states.size() ; k ++) {
		if (lprod[k] >= threshold || rprod[k] >= threshold) {
			state_indexes[hap].push_back(states[k]);
			state_lprobs[hap].push_back(lprod[k]);
			state_rprobs[hap].push_back(rprod[k]);
		}
	}
}
