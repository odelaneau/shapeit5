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

#include <phaser/phaser_header.h>

using namespace std;

void * hmmcompute_callback(void * ptr) {
	phaser * S = static_cast< phaser * >( ptr );
	int id_job, id_thread;

	pthread_mutex_lock(&S->mutex_workers);
	id_thread = S->i_threads ++;
	pthread_mutex_unlock(&S->mutex_workers);

	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= S->G.n_samples) vrb.progress("  * Processing", (id_job+1)*1.0/S->G.n_samples);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->G.n_samples) S->hmmcompute(id_job, id_thread);
		else pthread_exit(NULL);
	}
}

void phaser::hmmcompute(int id_job, int id_thread) {
	//Mapping storage events
	vector < vector < unsigned int > > cevents;
	G.mapUnphasedOntoScaffold(id_job, cevents);

	//Viterbi paths
	vector < int > path0, path1;

	//Forward-Backward-Viterbi passes for hap0
	thread_hmms[id_thread]->setup(2*id_job+0);
	thread_hmms[id_thread]->viterbi(path0);
	double pf0 = thread_hmms[id_thread]->forward();
	thread_hmms[id_thread]->backward(cevents, path0);


	//Forward-Backward-Viterbi passes for hap1
	thread_hmms[id_thread]->setup(2*id_job+1);
	thread_hmms[id_thread]->viterbi(path1);
	double pf1 = thread_hmms[id_thread]->forward();
	thread_hmms[id_thread]->backward(cevents, path1);

	//Phase remaining unphased using viterbi [singletons, etc ...]
	G.phaseCoalescentViterbi(id_job, path0, path1, M);
}

void phaser::phase() {
	//STEP1: haplotype selection
	vrb.title("PBWT pass");
	H.initialize(V,	options["pbwt-modulo"].as < double > (),
			options["pbwt-mdr"].as < double > (),
			options["pbwt-depth-common"].as < int > (),
			options["pbwt-depth-rare"].as < int > (),
			options["pbwt-mac"].as < int > ());
	H.scanIBD2(V);
	H.select(V, G);

	//STEP2: HMM computations
	vrb.title("HMM computations");
	thread_hmms = vector < hmm_scaffold * > (nthreads);
	for(int t = 0; t < nthreads ; t ++) thread_hmms[t] = new hmm_scaffold(V, G, H, M);
	if (nthreads > 1) {
		i_jobs = i_threads = 0;
		for (int t = 0 ; t < nthreads ; t++) pthread_create( &id_workers[t] , NULL, hmmcompute_callback, static_cast<void *>(this));
		for (int t = 0 ; t < nthreads ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_samples ; i ++) {
		hmmcompute(i, 0);
		vrb.progress("  * Processing", (i+1)*1.0/G.n_samples);
	}
	for(int t = 0; t < nthreads ; t ++) delete thread_hmms[t];
	vrb.bullet("Processing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	//STEP3: MERGE BACK ALL TOGETHER
	G.merge_by_transpose_I2V();

	//VERBOSE
	vrb.bullet("Solving summary:");
	G.nmiss_total = G.nmiss_imputation + G.nmiss_families + G.nmiss_monomorphic;
	G.nhets_total = G.nhets_families + G.nhets_imputation + G.nhets_coalescent;
	vrb.bullet2("#Hets=" + stb.str(G.nhets_total) + " / Phased by HMM=" + stb.str(G.nhets_imputation) + ", PED=" + stb.str(G.nhets_families) + ", SING=" + stb.str(G.nhets_coalescent));
	vrb.bullet2("#Miss=" + stb.str(G.nmiss_total) + " / Imputed by HMM=" + stb.str(G.nmiss_imputation) + ", PED=" + stb.str(G.nmiss_families) + ", MONO=" + stb.str(G.nmiss_monomorphic));

}
