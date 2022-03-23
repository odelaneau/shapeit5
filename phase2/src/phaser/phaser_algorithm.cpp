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
#include <phaser/phaser_header.h>

#include <models/hmm_scaffold/hmm_scaffold_header.h>
#include <models/gibbs_sampler/gibbs_sampler_header.h>


void * hmmcompute_callback(void * ptr) {
	phaser * S = static_cast< phaser * >( ptr );
	int id_job;
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= S->G.n_samples) vrb.progress("  * Processing", (id_job+1)*1.0/S->G.n_samples);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->G.n_samples) S->hmmcompute(id_job);
		else pthread_exit(NULL);
	}
}

void phaser::hmmcompute(int id_job) {
	vector < bool > cevents;
	vector < cprobs > cstates0, cstates1;
	hmm_scaffold HMM0(2*id_job+0, V, G, H, M);
	hmm_scaffold HMM1(2*id_job+1, V, G, H, M);
	HMM0.forward();
	HMM1.forward();
	G.mapUnphasedOntoScaffold(id_job, cevents);
	//unsigned int nstates0 = HMM0.backward(cevents, cstates0, options["compress-threshold"].as < double > ());
	//unsigned int nstates1 = HMM1.backward(cevents, cstates1, options["compress-threshold"].as < double > ());
	unsigned int nstates0 = HMM0.backward(cevents, cstates0);
	unsigned int nstates1 = HMM1.backward(cevents, cstates1);

	if (nthreads > 1) pthread_mutex_lock(&mutex_workers);
	Kstored.push(nstates0);
	Kstored.push(nstates1);
	for (int e0 = 0 ; e0 < cstates0.size() ; e0 ++) P.Pstates.push_back(cstates0[e0]);
	for (int e1 = 0 ; e1 < cstates1.size() ; e1 ++) P.Pstates.push_back(cstates1[e1]);
	if (nthreads > 1) pthread_mutex_unlock(&mutex_workers);
}


void phaser::phase() {
	tac.clock();
	i_jobs = 0;

	//STEP1: haplotype selection
	H.select();

	//STEP2: HMM computations
	vrb.title("HMM computations");
	if (nthreads > 1) {
		for (int t = 0 ; t < nthreads ; t++) pthread_create( &id_workers[t] , NULL, hmmcompute_callback, static_cast<void *>(this));
		for (int t = 0 ; t < nthreads ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_samples ; i ++) {
		hmmcompute(i);
		vrb.progress("  * Processing", (i+1)*1.0/G.n_samples);
	}
	vrb.bullet("Processing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("#storage_per_sample = " + stb.str(Kstored.mean(), 2) + " +/- " + stb.str(Kstored.sd(), 2));

	//STEP3: Big transpose
	vrb.bullet("Size of compressed probabilities = " + stb.str(P.sizeBytes() * 1.0f / 1e9, 4) + " Gb");
	P.transpose();
	P.mapping(H.n_scaffold_variants);

	//STEP4: MCMC computations
	vrb.title("Gibbs sampler computations");
	int errorRare = 0, errorCommon = 0, totalRare = 0, totalCommon = 0;
	gibbs_sampler GS (G.n_samples, 20, 10);
	for (int vs = 1 ; vs < V.sizeScaffold() ; vs ++) {
		for (int vt = V.vec_scaffold[vs-1]->idx_full + 1 ; vt < V.vec_scaffold[vs]->idx_full ; vt ++) {
			float weight = (V.vec_full[vt]->cm - V.vec_scaffold[vs-1]->cm) / (V.vec_scaffold[vs]->cm - V.vec_scaffold[vs-1]->cm);
			if (V.vec_full[vt]->type == VARTYPE_RARE) {
				GS.loadRare(G, H, P, V.vec_full[vt]->idx_rare, weight);
				GS.iterate(errorRare, totalRare);
				GS.pushRare(G, V.vec_full[vt]->idx_rare);
			} else {
				assert(V.vec_full[vt]->type == VARTYPE_COMM);
				GS.loadCommon(G, H, P, V.vec_full[vt]->idx_common, weight);
				GS.iterate(errorCommon, totalCommon);
				GS.pushCommon(G, V.vec_full[vt]->idx_common);
			}
		}
		vrb.progress("  * Processing", (vs+1)*1.0/V.sizeScaffold());
	}
	vrb.bullet("Processing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("Error rate at rare = " + stb.str(errorRare*100.0/totalRare, 4));
	vrb.bullet("Error rate at common = " + stb.str(errorCommon*100.0/totalCommon, 4));
}
