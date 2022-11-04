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

#ifndef _IBD2_TRACKS_H
#define _IBD2_TRACKS_H

#include <utils/otools.h>

struct track {
	int ind, from, to;

	track(int _ind, int _from, int _to) {
		ind = _ind;
		from = _from;
		to = _to;
	}

	bool operator<(const track & rhs) const {
		if (ind < rhs.ind) return true;
		if (ind > rhs.ind) return false;
		return (from < rhs.from);
	}

	bool overlap (const track & rhs) const {
		return ((ind==rhs.ind) && (rhs.to >= from) && (rhs.from <= to));
	}

	void merge (const track & rhs) {
		from = std::min(from, rhs.from);
		to = std::max(to, rhs.to);
	}
};

class ibd2_tracks {
public:

	std::vector < std::vector < track > > IBD2;


	ibd2_tracks ();
	~ibd2_tracks ();
	void clear();
	void initialize(int);

	bool noIBD2(int, int, int);
	void pushIBD2(int, std::vector < track > &);

	int collapse(std::vector < track > &);
	void collapse();
};

#endif
