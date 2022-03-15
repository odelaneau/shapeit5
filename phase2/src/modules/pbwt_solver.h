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
#ifndef _PBWT_SOLVER2_H
#define _PBWT_SOLVER2_H

#include <utils/otools.h>
#include <containers/haplotype_set.h>

class pbwt_solver {
private:
	bitmatrix & H;
	unsigned int n_site, n_main_hap, n_ref_hap, n_total_hap, n_total, n_resolved;
	vector < vector < int > > pbwt_clusters, pbwt_indexes, pbwt_divergences;
	vector < char > Guess;
	vector < bool > Het, Mis, Amb;
	vector < float > scoreBit;

public:
	pbwt_solver(haplotype_set &);
	~pbwt_solver();
	void free();

	void sweep(genotype_set &);
};

#endif
