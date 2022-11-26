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

#ifndef _HMM_PARAMETERS_H
#define _HMM_PARAMETERS_H

#include <utils/otools.h>
#include <containers/variant_map.h>

class hmm_parameters {
public :
	//DATA
	unsigned int Neff, Nhap;
	std::vector < float > t, nt, cm;
	float ee, ed;

	//CONSTRUCTOR/DESTRUCTOR
	hmm_parameters();
	~hmm_parameters();

	//METHODS
	void initialise(variant_map &, unsigned int, unsigned int);
	float getTransProb(float dist_cm);
};

inline
float hmm_parameters::getTransProb(float dist_cm) {
	if (dist_cm <= 1e-7) return -1.0f * expm1f(-0.04 * Neff * 1e-7 / Nhap);
	else return -1.0f * expm1f(-0.04 * Neff * dist_cm / Nhap);
}

#endif
