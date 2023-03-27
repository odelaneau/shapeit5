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

#ifndef _WINDOW_SET_H
#define _WINDOW_SET_H

#include <utils/otools.h>

#include <objects/genotype/genotype_header.h>

#include <containers/variant_map.h>

class window {
public:
	int start_locus;
	int start_segment;
	int start_ambiguous;
	int start_missing;
	int start_transition;
	int stop_locus;
	int stop_segment;
	int stop_ambiguous;
	int stop_missing;
	int stop_transition;

	window() {
		start_locus = 0;
		start_segment = 0;
		start_ambiguous = 0;
		start_missing = 0;
		start_transition = 0;
		stop_locus = 0;
		stop_segment = 0;
		stop_ambiguous = 0;
		stop_missing = 0;
		stop_transition = 0;
	}

	~window() {
	}

	std::string toString() {
		std::string str="";
		str += "L=[" + stb.str(start_locus) + "->" + stb.str(stop_locus) + "]";
		return str;
	}

	int lengthBP(variant_map & V) {
		return V.vec_pos[stop_locus]->bp - V.vec_pos[start_locus]->bp;
	}

	float lengthCM(variant_map & V) {
		return V.vec_pos[stop_locus]->cm - V.vec_pos[start_locus]->cm;
	}
};

class window_set {
public:
	//
	std::vector < window > W;

	//
	window_set();
	~window_set();
	void clear();

	//
	int size();
	bool split(double, int, int, std::vector < int > &, std::vector < int > &, std::vector < double > &, std::vector < double > &, std::vector < int > &);
	int build (variant_map &, genotype *, float);
};

#endif
