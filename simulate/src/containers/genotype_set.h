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

#ifndef _GENOTYPE_SET_H
#define _GENOTYPE_SET_H

#include <containers/haplotype_set.h>

class genotype_set {
public:
	haplotype_set H;

	//MASK
	std::vector < std::vector < int > > IDX;
	std::vector < std::string > SID;
	std::vector < bool > DIP;

	//Output

	std::vector < std::vector < bool > > GEN;
	std::vector < std::vector < bool > > MIS;

	genotype_set(haplotype_set &);
	~genotype_set();

	//FAMILY
	void includeFamily(int n_twins, int n_obs_fam, int n_obs_fam_off, int n_hid_fam, int n_hid_fam_off);
	void fillFamily();
	void missingFamily(float);

	//HAPLOID
	void includeHaploid(int n_haploid, int n_diploid);
	void fillHaploid();
	void missingHaploid(float, float);

	//IO
	void writeValidation(std::string);
	void writeEstimation(std::string);
};

#endif
