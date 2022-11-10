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

#define TYPE_SNP	0
#define TYPE_INDEL	1

class haplotype_set {
public:

	//Counts
	unsigned int n_variants;

	//Validation Data
	std::vector < std::vector < bool > > Htrue;
	std::vector < std::vector < bool > > Missing;
	std::vector < std::vector < bool > > Phased;

	//Estimated Data
	std::vector < std::vector < bool > > Hesti;
	std::vector < std::vector < bool > > Hprob;
	std::vector < std::vector < bool > > Estimated;
	std::map < std::string, float > Vprob;
	std::vector < int > IDXesti;

	//Variant Data [Columns]
	std::vector < int > MAC;
	std::vector < bool > MinorAlleles;
	std::vector < int > Positions;
	std::vector < std::string > RSIDs;
	std::vector < std::string > REFs, ALTs;

	//Family Data [Rows]
	std::vector < std::string > vecSamples;
	std::map < std::string, int > mapSamples;
	std::vector < int > Mothers, Fathers;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();
	void clear();

	void push(std::string &);
	void readPedigrees(std::string, bool);
	void assumePhased();


};

#endif
