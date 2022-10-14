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

#ifndef _HAPLOTYPE_CHECKER_H
#define _HAPLOTYPE_CHECKER_H

#include <utils/otools.h>
#include <containers/haplotype_set.h>

class haplotype_checker {
public:
	//DATA
	haplotype_set & H;
	vector < vector < bool > > Errors;
	vector < vector < bool > > Checked;
	vector < vector < float > > Calib;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_checker(haplotype_set &, int);
	~haplotype_checker();

	//Routines
	void check();
	bool isSNP(string &, string &);

	//Summarize
	void writePerSample(string);
	void writePerVariant(string);
	void writePerFrequency(string);
	void writePerType(string);
	void writeBlock(string);
	void writeFlipSwitchErrorPerSample(string);
	void writeCalibration(string);
};

inline
bool haplotype_checker::isSNP(string & ref, string &alt) {
	bool bref = (ref == "A") || (ref == "T") || (ref == "G") || (ref == "C");
	bool balt = (alt == "A") || (alt == "T") || (alt == "G") || (alt == "C");
	return bref && balt;
}

#endif
