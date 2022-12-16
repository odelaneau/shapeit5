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
#include <containers/genotype_set/genotype_set_header.h>

class genotype_reader {
public:
	//DATA
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;

	//COUNTS
	unsigned int n_total_variants;
	unsigned int n_scaffold_variants;
	unsigned int n_rare_variants;
	unsigned int n_samples;
	std::vector < unsigned long > n_scaffold_genotypes;
	std::vector < unsigned long > n_rare_genotypes;

	//PARAMETERS
	std::string funphased;
	std::string fphased;
	std::string fbinary;
	std::string foutput;//for reporting
	std::string scaffold_region;
	int input_start;
	int input_stop;
	int nthreads;
	float minmaf;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_reader(haplotype_set &, genotype_set &, variant_map &);
	~genotype_reader();

	//PARAMS
	void setFilenames(std::string, std::string, std::string, std::string);
	void setThreads(int);
	void setRegions(std::string, int, int);
	void setMAF(float);

	//IO
	void scanGenotypesPlain();
	void readGenotypesPlain();
	void scanGenotypesSparse();
	void readGenotypesSparse();
	void allocateGenotypes();
};


#endif
