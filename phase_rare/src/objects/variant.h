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

#ifndef _VARIANT_H
#define _VARIANT_H

#include <utils/otools.h>

#define VARTYPE_SCAF	0
#define VARTYPE_COMM	1
#define VARTYPE_RARE	2

class variant {
public :
	//DATA
	std::string chr;
	int bp;
	std::string id;
	std::string ref;
	std::string alt;
	bool minor;
	double cm;
	unsigned int cref;
	unsigned int calt;
	unsigned int cmis;
	//NEW DATA FIELDS
	char type;
	int idx_full;
	int idx_scaffold;
	int idx_common;
	int idx_rare;


	//CONSTRUCTOR/DESTRUCTOR
	variant(std::string & chr, int bp, std::string & id, std::string & ref, std::string & alt, bool, char);
	~variant();

	bool isSingleton();
	bool isMonomorphic();
	bool isSNP();
	unsigned int getMAC();
	double getMDR();
	double getMAF();
	double getAF();
};

#endif
