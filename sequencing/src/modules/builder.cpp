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

#include <modules/builder.h>

builder::builder(genotype_set & _G, int n_thread): G(_G) {
	this->n_thread = n_thread;
	if (n_thread > 1) {
		i_workers = 0;
		id_workers = vector < pthread_t > (n_thread);
		pthread_mutex_init(&mutex_workers, NULL);
	}
}

builder::~builder() {
	if (n_thread > 1) {
		i_workers = 0;
		pthread_mutex_destroy(&mutex_workers);
		id_workers.clear();
	}
}

void * builder_callback(void * ptr) {
	builder * B = static_cast< builder * >( ptr );
	for(;;) {
		pthread_mutex_lock( &B->mutex_workers );
		int curr_ind_to_process = B->i_workers++;
		pthread_mutex_unlock( &B->mutex_workers);
		if (curr_ind_to_process < B->G.n_target_samples) B->build(curr_ind_to_process);
		else pthread_exit(NULL);
	}
	return NULL;
}

void builder::build(int ind) {
	G.vecG[ind]->build();
}

void builder::build() {
	tac.clock();
	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, builder_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i  <  G.n_target_samples ; i ++) build(i);
	long int n_segments = G.numberOfSegments();
	vrb.bullet("Build genotype graphs [seg=" + stb.str(n_segments) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)");
}

