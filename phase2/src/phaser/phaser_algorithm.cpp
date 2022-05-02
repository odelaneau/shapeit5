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

#include <models/gibbs_sampler/gibbs_sampler_header.h>
#include <models/pileup_caller/pileup_caller_header.h>

#define ALLOC_CHUNK 40000000

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
	vector < cstate > cstates0, cstates1;

	//Mapping storage events
	vector < bool > cevents;
	G.mapUnphasedOntoScaffold(id_job, cevents);

	//HMM compute
	thread_hmms[id_thread]->setup(2*id_job+0);
	thread_hmms[id_thread]->forward();
	thread_hmms[id_thread]->backward(cevents, cstates0);
	thread_hmms[id_thread]->setup(2*id_job+1);
	thread_hmms[id_thread]->forward();
	thread_hmms[id_thread]->backward(cevents, cstates1);

	//Storage of compressed probabilities [Mutex protected]
	if (nthreads > 1) pthread_mutex_lock(&mutex_workers);

	unsigned long int allocated_size = P.Pstates.capacity();
	unsigned long int necessary_size = P.Pstates.size()+cstates0.size()+cstates1.size();
	unsigned long int requested_size = 0;

	if (necessary_size > allocated_size) {
		if (id_job < (G.n_samples/5)) {
			requested_size = allocated_size + ALLOC_CHUNK;
		} else {
			requested_size = allocated_size + (G.n_samples - id_job + 1) * statCS.mean() * 1.1f;
		}
		P.Pstates.reserve(requested_size);
	}

	for (int e = 0 ; e < cstates0.size() ; e ++) P.Pstates.push_back(cstates0[e]);
	for (int e = 0 ; e < cstates1.size() ; e ++) P.Pstates.push_back(cstates1[e]);
	statCS.push(cstates0.size()+cstates1.size());

	if (nthreads > 1) pthread_mutex_unlock(&mutex_workers);
}


void * gibbscompute_callback(void * ptr) {
	phaser * S = static_cast< phaser * >( ptr );
	int id_job;
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->nthreads) S->gibbscompute(id_job);
		else pthread_exit(NULL);
	}
}

void phaser::gibbscompute(int id_job) {
	gibbs_sampler GS (G.n_samples, options["mcmc-iterations"].as < int > (), options["mcmc-burnin"].as < int > ());
	for (int v = 0 ; v < thread_data[id_job].size() ; v ++) {
		unsigned int _n_common_yphased = 0, _n_common_nphased = 0, _n_rare_yphased = 0, _n_rare_nphased = 0;
		int vt = thread_data[id_job][v].first;
		float weight = thread_data[id_job][v].second;
		if (V.vec_full[vt]->type == VARTYPE_RARE) {
			GS.loadRare(G, H, P, V.vec_full[vt]->idx_rare, weight);
			V.vec_full[vt]->isSNP()?GS.setHQ():GS.setLQ();
			GS.iterate();
			GS.pushRare(G, V.vec_full[vt]->idx_rare, _n_rare_yphased, _n_rare_nphased, 0.8f);
		} else {
			assert(V.vec_full[vt]->type == VARTYPE_COMM);
			GS.loadCommon(G, H, P, V.vec_full[vt]->idx_common, weight);
			V.vec_full[vt]->isSNP()?GS.setHQ():GS.setLQ();
			GS.iterate();
			GS.pushCommon(G, V.vec_full[vt]->idx_common, _n_common_yphased, _n_common_nphased, 0.249f);
		}

		if (nthreads > 1) pthread_mutex_lock(&mutex_workers);
		doneSite++;
		n_common_yphased += _n_common_yphased;
		n_common_nphased +=_n_common_nphased;
		n_rare_yphased += _n_rare_yphased;
		n_rare_nphased += _n_rare_nphased;
		vrb.progress("  * Processing", doneSite*1.0/totalSite);
		if (nthreads > 1) pthread_mutex_unlock(&mutex_workers);
	}
}

void phaser::phase() {
	//STEP0: SINGLETONS AND MONOMORPHIC
	vrb.title("Pre-processing");
	G.imputeMonomorphic();
	//G.randomizeSingleton();

	//STEP1: haplotype selection
	vrb.title("PBWT pass");
	H.initialize(V,	options["pbwt-modulo"].as < double > (),
			options["pbwt-mdr"].as < double > (),
			options["pbwt-depth-common"].as < int > (),
			options["pbwt-depth-rare"].as < int > (),
			options["pbwt-mac"].as < int > ());

	G.transpose();
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
	vrb.bullet("#states = " + stb.str(P.Pstates.size()) + " / #reserved = " + stb.str(P.Pstates.capacity()) + " / Memory = " + stb.str(P.Pstates.capacity() * sizeof(cstate)/1e9, 2) + "Gg");

	//STEP3: Big transpose
	P.transpose();
	P.mapping(H.n_scaffold_variants);

	//STEP4: MCMC loading
	vector < pair < int, float > > thread_data_serialized;
	totalSite = 0; doneSite = 0;
	for (int vt = 0 ; vt < V.vec_scaffold[0]->idx_full ; vt ++) {
		thread_data_serialized.push_back(pair < int, float > (vt, 0.5f));
		totalSite++;
	}
	for (int vs = 1 ; vs < V.sizeScaffold() ; vs ++) {
		for (int vt = V.vec_scaffold[vs-1]->idx_full + 1 ; vt < V.vec_scaffold[vs]->idx_full ; vt ++) {
			float weight = 0.5f;
			if ((V.vec_scaffold[vs]->cm - V.vec_scaffold[vs-1]->cm) > 1e-7) {
				weight = (V.vec_full[vt]->cm - V.vec_scaffold[vs-1]->cm) / (V.vec_scaffold[vs]->cm - V.vec_scaffold[vs-1]->cm);
				weight = max (0.00001f, weight);
				weight = min (0.99999f, weight);
			}
			thread_data_serialized.push_back(pair < int, float > (vt, weight));
			totalSite++;
		}
	}
	for (int vt = V.vec_scaffold.back()->idx_full + 1 ; vt < V.sizeFull() ; vt ++) {
		thread_data_serialized.push_back(pair < int, float > (vt, 0.5f));
		totalSite++;
	}

	//STEP5: MCMC computations
	vrb.title("Gibbs sampler computations");
	n_common_yphased = 0;
	n_common_nphased = 0;
	n_rare_yphased = 0;
	n_rare_nphased = 0;
	thread_data = vector < vector < pair < int, float > > > (nthreads);
	for(int j = 0; j < thread_data_serialized.size() ; j ++) thread_data[j%nthreads].push_back(thread_data_serialized[j]);
	if (nthreads > 1) {
		i_jobs = 0;
		for (int t = 0 ; t < nthreads ; t++) pthread_create( &id_workers[t] , NULL, gibbscompute_callback, static_cast<void *>(this));
		for (int t = 0 ; t < nthreads ; t++) pthread_join( id_workers[t] , NULL);
	} else gibbscompute(0);
	vrb.bullet("Processing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	if (G.n_common_variants) vrb.bullet2("#common_phased = " + stb.str(n_common_yphased) + " / " + stb.str(n_common_yphased+n_common_nphased) + " (" + stb.str(n_common_yphased * 100.0 / (n_common_yphased+n_common_nphased), 2) + "%)");
	vrb.bullet2("#rare_phased = " + stb.str(n_rare_yphased) + " / " + stb.str(n_rare_yphased+n_rare_nphased) + " (" + stb.str(n_rare_yphased * 100.0 / (n_rare_yphased+n_rare_nphased), 2) + "%)");

	//STEP6: PHASE REMAINING HETS
	vrb.title("Solving remaining hets using PBWT");
	H.solve(V, G);


}
