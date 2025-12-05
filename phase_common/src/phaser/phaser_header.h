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

#include <models/haplotype_segment_single.h>
#include <models/haplotype_segment_double.h>

#include <containers/genotype_set.h>
#include <containers/conditioning_set/conditioning_set_header.h>
#include <containers/variant_map.h>

#define STAGE_BURN	0
#define STAGE_PRUN	1
#define STAGE_MAIN	2

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

	//PBWT
	bool pbwt_auto;
	int pbwt_depth;
	double pbwt_modulo;

	//MULTI-THREADING
	int i_workers, i_jobs;
	std::vector < pthread_t > id_workers;
	pthread_mutex_t mutex_workers;
	std::vector < compute_job > threadData;

	//MCMC
	std::vector < unsigned int > iteration_types;
	std::vector < unsigned int > iteration_counts;
	unsigned int iteration_stage;
	int n_underflow_recovered_summing;
	int n_underflow_recovered_precision;

	//STATS
	stats1D statH,statS;
	std::vector < double > storedKsizes;

	//CONSTRUCTOR
	phaser();
	~phaser();

	//METHODS
	void phase();
	void phaseWindow(int, int);
	void phaseWindow();

	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void parse_iteration_scheme(std::string);
	std::string get_iteration_scheme();
	void check_options();
	void verbose_options();
	void verbose_files();

	//
	void read_files_and_initialise();
	void phase(std::vector < std::string > &);
	void write_files_and_finalise();
};


#endif


