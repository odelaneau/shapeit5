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

#ifndef _GENOTYPE_SET_H
#define _GENOTYPE_SET_H

#include <utils/otools.h>

#include <objects/genotype/genotype_header.h>
#include <containers/variant_map.h>
#include <io/pedigree_reader.h>


class genotype_set {
public:
	//DATA
	int n_site, n_ind;							//Number of variants, number of individuals
	std::vector < genotype * > vecG;					//Vector of genotype graphs
	std::vector < genotype * > vecFathers;			//Points to fathers, NULL otherwise
	std::vector < genotype * > vecMothers;			//Points to mothers, NULL otherwise

	//CONSTRUCTOR/DESTRUCTOR
	genotype_set();
	~genotype_set();
	void allocate(unsigned long, unsigned long);

	//METHODS
	void imputeMonomorphic(variant_map &);		//Impute to REF monomorphic variants
	unsigned int largestNumberOfTransitions();	//Get the number of transitions in the larger genotype graph. Used to initialize memory space for multi-threading.
	unsigned int largestNumberOfMissings();		//Get the number of transitions in the larger genotype graph. Used to initialize memory space for multi-threading.
	unsigned long numberOfSegments();			//Total number of segments across all genotype graphs (used for verbose).
	void solve();								//
	void scaffoldUsingPedigrees(pedigree_reader &);
	void resetHaploidHeterozgotes(std::vector < std::string > &);

};

#endif
