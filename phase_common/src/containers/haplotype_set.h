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

#ifndef _HAPLOTYPE_SET_H
#define _HAPLOTYPE_SET_H

#include <utils/otools.h>

#include <containers/bitmatrix.h>
#include <containers/genotype_set.h>
#include <containers/variant_map.h>

#include <io/pedigree_reader.h>


class haplotype_set {
public:
	//Haplotype Data
	bitmatrix H_opt_hap;		// Bit matrix of haplotypes (haplotype first). Transposed version of H_opt_var.
	bitmatrix H_opt_var;		// Bit matrix of haplotypes (variant first). Transposed version of H_opt_hap.
	unsigned long n_site;		// #variants
	unsigned long n_hap;		// #haplotypes
	unsigned long n_ind;		// #individuals

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();
	void clear();
	void allocate(unsigned long, unsigned long, unsigned long);

	//Haplotype routines
	void updateHaplotypes(genotype_set & G, bool first_time = false);
	void transposeHaplotypes_H2V(bool full, bool verbose = true);
	void transposeHaplotypes_V2H(bool full, bool verbose = true);
};

#endif
