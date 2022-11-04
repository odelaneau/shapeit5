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

#ifndef _MENDEL_SOLVER_H
#define _MENDEL_SOLVER_H

#include <utils/otools.h>
#include <containers/haplotype_set.h>

class mendel_solver {
public:

	//DATA
	haplotype_set & H;
	std::vector < std::vector < bool > > Errors;
	std::vector < int > CountsD0, CountsD1, CountsD2, CountsT00, CountsT01, CountsT02, CountsT10, CountsT11, CountsT12;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	mendel_solver(haplotype_set &);
	~mendel_solver();

	//Routines
	void solveT(int locus, int cidx, int fidx, int midx);
	void solveD(int locus, int cidx, int pidx, bool father, bool singleton);
	void solve(bool);
	void countT(int locus, int cidx, int fidx, int midx);
	void countD(int locus, int cidx, int pidx);
	void count();
	void set();

	//Summarize
	void writePerSample(std::string);
	void writePerVariant(std::string);
	void writeImbalance(std::string);
	void writePedigree(std::string);
};

#endif
