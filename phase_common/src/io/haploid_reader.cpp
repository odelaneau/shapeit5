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

#include <io/haploid_reader.h>

using namespace std;

haploid_reader::haploid_reader() {
}

haploid_reader::~haploid_reader() {
	vector < string > ().swap(samples);
}

void haploid_reader::readHaploidFile(string fmales) {
	tac.clock();
	string buffer;
	input_file fd_males(fmales);
	if (fd_males.fail()) vrb.error("Cannot open haploid [e.g. chrX males] file");
	while (getline(fd_males, buffer)) samples.push_back(buffer);
	fd_males.close();
	vrb.bullet("Haploid parsing [n=" + stb.str(samples.size()) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
