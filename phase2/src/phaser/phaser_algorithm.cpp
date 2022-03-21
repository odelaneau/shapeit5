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
	hmm_scaffold HMM0(2*id_job+0, V, G, H, M);
	hmm_scaffold HMM1(2*id_job+1, V, G, H, M);
	HMM0.forward();
	HMM1.forward();
	G.mapUnphasedOntoScaffold(id_job, cevents);
	HMM0.backward(cevents, P.Phap_states[2*id_job+0], P.Phap_indexes[2*id_job+0], options["compress-threshold"].as < double > ());
	HMM1.backward(cevents, P.Phap_states[2*id_job+1], P.Phap_indexes[2*id_job+1], options["compress-threshold"].as < double > ());

	pthread_mutex_lock(&mutex_workers);
	Kstored.push(P.Phap_states[2*id_job+0].size());
	Kstored.push(P.Phap_states[2*id_job+1].size());
	pthread_mutex_unlock(&mutex_workers);
}


void phaser::phase() {
	tac.clock();
	int n_threads = options["thread"].as < int > ();
	i_jobs = 0;

	//STEP1: haplotype selection
	H.select();

	//STEP2: HMM computations
	vrb.title("HMM computations");
	if (n_threads > 1) {
		for (int t = 0 ; t < n_threads ; t++) pthread_create( &id_workers[t] , NULL, hmmcompute_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_threads ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_samples ; i ++) {
		hmmcompute(i);
		vrb.progress("  * Processing", (i+1)*1.0/G.n_samples);
	}
	vrb.bullet("Processing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("#storage_per_sample=" + stb.str(Kstored.mean(), 2));

	//STEP3: MCMC computations
	//... TO DO ...
}
