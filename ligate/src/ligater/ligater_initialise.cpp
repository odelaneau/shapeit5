/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
 * Copyright (C) 2022-2023 Simone Rubinacci
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

#include <ligater/ligater_header.h>

using namespace std;

void ligater::read_files_and_initialise() {

	//step0: Initialize seed & other
	rng.setSeed(options["seed"].as < int > ());

	//step1:read filenames
	string buffer;
	string filelist = options["input"].as < string > ();
	vrb.title("Read filenames in [" + filelist + "]");
	input_file fd(filelist);
	while (getline(fd, buffer)) filenames.push_back(buffer);
	vrb.bullet("#files = " + stb.str(filenames.size()));
	if (filenames.size() == 0) vrb.error("No filenames in input file.");
	nfiles = filenames.size();

	//step2: read PED file
	if (options.count("pedigree")) {
		vrb.title("Read filenames in [" + filelist + "]");
		string buffer;
		vector < string > tokens;
		input_file fd_ped(options["pedigree"].as < string > ());
		if (fd_ped.fail()) vrb.error("Cannot open pedigree file");
		while (getline(fd_ped, buffer)) {
			stb.split(buffer, tokens);
			if (tokens.size() < 3) vrb.error("Problem in pedigree file; line should have at least 3 columns");
			kids.push_back(tokens[0]);
			fathers.push_back(tokens[1]);
			mothers.push_back(tokens[2]);
		}
		fd_ped.close();
		vrb.bullet("PED file parsing [n=" + stb.str(kids.size()) + "]");
	}
}
