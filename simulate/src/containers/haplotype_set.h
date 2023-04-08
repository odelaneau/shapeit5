/*******************************************************************************
 * Copyright (C) 2020 Olivier Delaneau, University of Lausanne
 * Copyright (C) 2020 Simone Rubinacci, University of Lausanne
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

class haplotype_set {
public:

	//SAMPLES
	std::vector < std::vector < bool > > HAP;
	std::vector < std::string > SID;
	unsigned int n_hap;

	//VAR
	std::set < int > filter_positions;
	unsigned int n_var;
	std::vector < std::string > VID;
	std::vector < int > POS;
	std::vector < std::string > REF;
	std::vector < std::string > ALT;
	std::vector < std::string > CHR;
	std::vector < float > MAF;

	//CONSTRUCTORS/DESCTRUCTORS
	haplotype_set();
	~haplotype_set();

	//IO
	void readPositions(std::string fpos);
	void readHaplotypes(std::string fdata, std::string region, bool haploid);

	//ROUTINES
	bool checkPos(int);
	void computeMAF();
};

#endif
