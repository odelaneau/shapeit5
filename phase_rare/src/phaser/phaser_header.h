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

#ifndef _PHASER_H
#define _PHASER_H

#include <utils/otools.h>
#include <objects/hmm_parameters.h>

#include <containers/genotype_set/genotype_set_header.h>
#include <containers/state_set.h>
#include <containers/conditioning_set/conditioning_set_header.h>
#include <containers/variant_map.h>

#include <models/hmm_scaffold/hmm_scaffold_header.h>


class phaser {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//INTERNAL DATA
	conditioning_set H;
	genotype_set G;
	hmm_parameters M;
	variant_map V;
	state_set P;

	//MULTI-THREADING
	int i_jobs, i_threads, nthreads;
	std::vector < pthread_t > id_workers;
	pthread_mutex_t mutex_workers;
	std::vector < std::vector < std::pair < int, float > > > thread_data;
	std::vector < hmm_scaffold * > thread_hmms;

	//STATS
	int totalSite, doneSite;
	unsigned long int n_common_yphased;
	unsigned long int n_common_nphased;
	unsigned long int n_rare_yphased;
	unsigned long int n_rare_nphased;
	stats1D statCS;

	//GENOMIC REGION
	std::string chrid;
	int input_start;
	int input_stop;
	int scaffold_start;
	int scaffold_stop;
	std::string input_gregion;
	std::string scaffold_gregion;

	//CONSTRUCTOR
	phaser();
	~phaser();

	//METHODS
	void hmmcompute(int, int);
	void phase();


	//PARAMETERS
	void buildCoordinates();
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//
	void read_files_and_initialise();
	void phase(std::vector < std::string > &);
	void write_files_and_finalise();
};


#endif


