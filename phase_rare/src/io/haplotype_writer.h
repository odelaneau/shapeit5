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

#ifndef _HAPLOTYPE_WRITER_H
#define _HAPLOTYPE_WRITER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>
#include <containers/genotype_set/genotype_set_header.h>


class haplotype_writer {
public:
	//DATA
	int nthreads;
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;
	int input_start;
	int input_stop;

	//CONSTRUCTORS/DESCTRUCTORS
	haplotype_writer(haplotype_set &, genotype_set &, variant_map &, int);
	~haplotype_writer();
	void setRegions(int _input_start, int _input_stop);

	//IO
	void writeHaplotypesVCF(std::string foutput, std::string ifile, bool addPP);
	void writeHaplotypesXCF(std::string foutput, std::string ifile, std::string format);
};

#endif
