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

#ifndef _CONDITIONING_SET_H
#define _CONDITIONING_SET_H

#include <utils/otools.h>

#include <containers/haplotype_set.h>
#include <containers/ibd2_tracks.h>

class conditioning_set : public haplotype_set {
public:

	//VARIANT INDEXING
	std::vector < int > sites_pbwt_mthreading;
	std::vector < bool > sites_pbwt_evaluation;
	std::vector < bool > sites_pbwt_selection;
	std::vector < int > sites_pbwt_grouping;
	//std::vector < int > sites_pbwt_storage;
	std::vector < int > starts_pbwt_mthreading;
	unsigned int sites_pbwt_ngroups;

	//PARAMETERS FOR PBWT
	int depth, nthread;

	//IBD2 TRACKS
	ibd2_tracks Kbanned;

	//STATE DATA
	std::vector < int > indexes_pbwt_neighbour;

	//SOLVER DATA
	std::vector < float > scoreBit;

	//MULTI-THREADING
	int i_worker, i_job, d_job;
	pthread_mutex_t mutex_workers;
	std::vector < pthread_t > id_workers;

	//CONSTRUCTOR/DESTRUCTOR
	conditioning_set();
	~conditioning_set();
	void initialize(variant_map & V, float _modulo_selection, float _modulo_multithreading, float _mdr, int _depth, int _mac, int _nthread);

	//VARIANT PROCESSING
	bool split(variant_map & V, float min_length, int left_index, int right_index, std::vector < int > & output);

	//STATES PROCESSING
	void store(int l, std::vector < int > & A, std::vector < int > & C);
	void select(int chunk);
	void select();
	void transposePBWTneighbours();

	//PBWT PHASING SWEEP
	void solve(int chunk, genotype_set *);
	void solve(genotype_set *);

};

#endif
