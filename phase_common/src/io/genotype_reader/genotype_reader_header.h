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

	//COUNTS
	unsigned long n_variants;
	unsigned long n_main_samples;
	unsigned long n_ref_samples;
	std::vector < unsigned long > n_genotypes;

	//PARAMETERS
	std::vector< std::string > filenames;
	std::vector< char > panels;
	std::string region;
	int threads;
	float filter_min_maf;
	bool filter_snp_only;
	std::vector < bool > variant_mask;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_reader(haplotype_set &, genotype_set &, variant_map &);
	~genotype_reader();

	//PARAMS
	void setFilterMAF(float);
	void setFilterSNP();
	void setMainFilename(std::string);
	void addReferenceFilename(std::string);
	void addScaffoldFilename(std::string);
	void setThreads(int);
	void setRegion(std::string);

	//IO
	void scanGenotypes();
	void readGenotypes();
	void allocateGenotypes();
};


#endif
