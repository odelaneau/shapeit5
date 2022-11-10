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

#include <io/pedigree_reader.h>

using namespace std;

pedigree_reader::pedigree_reader() {
}

pedigree_reader::~pedigree_reader() {
	vector < string > ().swap(fathers);
	vector < string > ().swap(mothers);
	vector < string > ().swap(kids);
}

void pedigree_reader::readPedigreeFile(string fped) {
	tac.clock();
	string buffer;
	vector < string > tokens;
	int line = 0;
	input_file fd_ped(fped);
	if (fd_ped.fail()) vrb.error("Cannot open pedigree file");
	while (getline(fd_ped, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() != 3) vrb.error("Problem in pedigree file; each line should have 3 columns");
		kids.push_back(tokens[0]);
		fathers.push_back(tokens[1]);
		mothers.push_back(tokens[2]);
		line ++;
	}
	fd_ped.close();
	vrb.bullet("PED parsing [n=" + stb.str(line) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
