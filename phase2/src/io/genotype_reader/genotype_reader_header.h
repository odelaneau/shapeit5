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
	variant_map & V;

	//COUNTS
	unsigned int n_scaffold_variants;
	unsigned int n_rare_variants;
	unsigned int n_common_variants;
	unsigned int n_samples;
	vector < unsigned long > n_scaffold_genotypes;
	vector < unsigned long > n_common_genotypes;
	vector < unsigned long > n_rare_genotypes;

	//PARAMETERS
	string funphased;
	string fphased;
	string region;
	int threads;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_reader(haplotype_set &, variant_map &);
	~genotype_reader();

	//PARAMS
	void setFilenames(string, string);
	void setThreads(int);
	void setRegion(string);

	//IO
	void scanGenotypes();
	void readGenotypes();
	void allocateGenotypes();
};


#endif
