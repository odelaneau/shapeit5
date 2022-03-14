/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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
#ifndef _GENOTYPE_READER_H
#define _GENOTYPE_READER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>

class genotype_reader {
public:
	//DATA
	int nthreads;
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;
	string region;

	//PARAMS
	bool use_PS_field;
	float maf_threshold;

	//COUNTS
	unsigned int n_common_variants;
	unsigned int n_rare_variants;
	unsigned int n_total_variants;
	unsigned int n_target_samples;
	unsigned int n_reference_samples;
	unsigned int n_total_samples;
	vector < unsigned long > n_genotypes;

	//FLAGS
	vector < bool > flag_common;

	//PHASESETS
	unordered_map < int, int > PSmap;
	vector < int > PScodes;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_reader(haplotype_set &, genotype_set &, variant_map &, string , bool _use_PS_field, int _nthreads);
	~genotype_reader();

	//IO
	void scanGenotypes(string funphased);
	void scanGenotypes(string funphased, string fphased);
	void allocateGenotypes();
	void readGenotypes0(string);
	void readGenotypes1(string, string);
	//void readGenotypes2(string, string);
	//void readGenotypes3(string, string, string);
	void setPScodes(int * ps_arr, int nps);
};

#endif
