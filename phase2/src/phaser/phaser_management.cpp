////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <phaser/phaser_header.h>

phaser::phaser() {
}

phaser::~phaser() {
	id_workers.clear();
	threadData.clear();
	iteration_types.clear();
	iteration_counts.clear();
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

void phaser::parse_iteration_scheme(string str_iter) {
	vector < string > tokens;
	if (stb.split(str_iter, tokens, ",") < 1) vrb.error("Impossible to parse iteration scheme [" + str_iter + "]");
	iteration_types.clear();
	iteration_counts.clear();
	for (int stage = 0 ; stage < tokens.size() ; stage ++) {
		int n_iter = std::stoi(tokens[stage]);
		string t_iter = tokens[stage].substr(tokens[stage].size()-1);
		if (n_iter < 0) vrb.error("Incorrect number of iterations [" + stb.str(n_iter) + "]");
		iteration_counts.push_back(n_iter);
		if (t_iter == "b") iteration_types.push_back(STAGE_BURN);
		else if (t_iter == "p") iteration_types.push_back(STAGE_PRUN);
		else if (t_iter == "m") iteration_types.push_back(STAGE_MAIN);
		else vrb.error("Unrecognized iteration type [" + t_iter + "]");
	}
}

string phaser::get_iteration_scheme() {
	int n_total_iter = 0;
	for (int s = 0 ; s < iteration_counts.size() ; s ++) n_total_iter += iteration_counts[s];

	string str = stb.str(n_total_iter) + " iterations [";
	for (int s = 0 ; s < iteration_counts.size() ; s ++) {
		str += (s?" + ":"") + stb.str(iteration_counts[s]);
		switch (iteration_types[s]) {
		case STAGE_BURN:	str+= "b"; break;
		case STAGE_PRUN:	str+= "p"; break;
		case STAGE_MAIN:	str+= "m"; break;
		}
	}
	str += "]";
	return str;
}

