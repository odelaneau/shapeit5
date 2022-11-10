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

#include <phaser/phaser_header.h>

using namespace std;

phaser::phaser() {
}

phaser::~phaser() {
	id_workers.clear();
}

void phaser::phase(vector < string > & args) {
	declare_options();
	parse_command_line(args);
	check_options();
	verbose_files();
	verbose_options();
	read_files_and_initialise();
	phase();
	write_files_and_finalise();
}

void phaser::buildCoordinates() {
        string scaffold_region = options["scaffold-region"].as < string > ();
        string input_region = options["input-region"].as < string > ();
        vrb.title("Parsing specified genomic regions");
        vector < string > scaffold_t1, scaffold_t2;
        vector < string > input_t1, input_t2;
        int scaffold_ret = stb.split(scaffold_region, scaffold_t1, ":");
        int input_ret = stb.split(input_region, input_t1, ":");
        if (scaffold_ret != 2) vrb.error("Input region needs to be specificied as chrX:Y-Z (chromosome ID cannot be extracted)");
        if (input_ret != 2) vrb.error("Scaffold region needs to be specificied as chrX:Y-Z (chromosome ID cannot be extracted)");
        chrid = scaffold_t1[0];
        if (chrid != input_t1[0]) vrb.error("Chromosome IDs in scaffold and input regions are different!");
        scaffold_ret = stb.split(scaffold_t1[1], scaffold_t2, "-");
        input_ret = stb.split(input_t1[1], input_t2, "-");
        if (scaffold_ret != 2) vrb.error("Input region needs to be specificied as chrX:Y-Z (genomic positions cannot be extracted)");
        if (input_ret != 2) vrb.error("input region needs to be specificied as chrX:Y-Z (genomic positions cannot be extracted)");
        scaffold_start = atoi(scaffold_t2[0].c_str());
        scaffold_stop = atoi(scaffold_t2[1].c_str());
        input_start = atoi(input_t2[0].c_str());
        input_stop = atoi(input_t2[1].c_str());
        if (scaffold_start >= scaffold_stop) vrb.error("Input genomic region coordinates are incorrect (start >= stop)");
        if (input_start >= input_stop) vrb.error("input genomic region coordinates are incorrect (start >= stop)");
        if (scaffold_start > input_start) vrb.error("Input/input genomic region coordinates are imcompatible (scaffold_start > input_start)");
        if (scaffold_stop < input_stop) vrb.error("Input/input genomic region coordinates are incompatible (scaffold_stop < input_stop)");
        if (scaffold_start < 0) vrb.error("Input genomic region coordinates are incorrect (scaffold_start < 0)");
        if (input_start < 0) vrb.error("Input genomic region coordinates are incorrect (input_start < 0)");
        scaffold_gregion = chrid + ":" + stb.str(scaffold_start) + "-" + stb.str(scaffold_stop);
        input_gregion = chrid + ":" + stb.str(input_start) + "-" + stb.str(input_stop);
        vrb.bullet("Scaffold region  [" + scaffold_gregion + "]");
        vrb.bullet("Input region  [" + input_gregion + "]");
}
