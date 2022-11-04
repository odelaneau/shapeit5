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

#define _GLOBAL
#include <objects/variant.h>

using std::string;

variant::variant(string & chr, int bp, string & id, string & ref, string & alt, bool _minor, char _type) {
	this->chr = chr;
	this->bp = bp;
	this->id = id;
	this->ref = ref;
	this->alt = alt;
	this->cref = 0;
	this->calt = 0;
	this->cmis = 0;
	type = _type ;
	minor = _minor;
	cm = -1;
	idx_full = -1;
	idx_scaffold = -1;
	idx_common = -1;
	idx_rare = -1;
}

variant::~variant() {
}

unsigned int variant::getMAC() {
	return std::min(cref, calt);
}

double variant::getMAF() {
	double freq = cref * 1.0f / (cref + calt);
	return std::min(freq, 1.0f-freq);
}

double variant::getAF() {
	return cref * 1.0f / (cref + calt);
}

double variant::getMDR() {
	return cmis * 1.0 / (cref + calt + cmis);
}

bool variant::isSingleton() {
	return (calt == 1 || cref == 1);
}

bool variant::isMonomorphic() {
	return (calt == 0 || cref == 0);
}

bool variant::isSNP() {
	bool bref = (ref == "A") || (ref == "T") || (ref == "G") || (ref == "C");
	bool balt = (alt == "A") || (alt == "T") || (alt == "G") || (alt == "C");
	return (bref && balt);
}


