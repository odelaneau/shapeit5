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

#include <switcher/switcher_header.h>

#include <io/haplotype_reader.h>

using std::string;

void switcher::read_files_and_initialise() {
	//step1: Read input files

	if (options.count("frequency"))
		haplotype_reader(H, options["region"].as < string > (), options["min-pp"].as < double > (), options["thread"].as < int > ()).readHaplotypes(options["validation"].as < string > (), options["estimation"].as < string > (), options["frequency"].as < string > (), options.count("dupid"));
	else
		haplotype_reader(H, options["region"].as < string > (), options["min-pp"].as < double > (), options["thread"].as < int > ()).readHaplotypes(options["validation"].as < string > (), options["estimation"].as < string > (), options.count("dupid"));

	//step2: read pedigrees if necessary
	if (options.count("pedigree")) H.readPedigrees(options["pedigree"].as < std::string > (), options.count("dupid"));
	else H.assumePhased();
}
