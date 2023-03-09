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

#include "../../versions/versions.h"

#include <concater/concater_header.h>

using namespace std;

void concater::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of thread used for VCF/BCF (de-)compression");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input-sparse-list", bpo::value< string >(), "Input genotype data in sparse VCF/BCF format")
			("region", bpo::value< string >(), "Region to be considered in --input");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output-sparse", bpo::value< string >(), "Output genotype data in sparse VCF/BCF format")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void concater::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[SHAPEIT5] Concat multiple sparse VCF/BCF");
	vrb.bullet("Authors       : Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 5." + string(SCFTLS_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());
}

void concater::check_options() {
	int ninput = options.count("input-sparse-list");
	int noutput = options.count("output-sparse");

	if (ninput != 1) vrb.error("Use --input-sparse-list as input.");
	if (noutput != 1) vrb.error("Use --output-sparse as output.");

	if (!options.count("region"))
		vrb.error("--region missing");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");
}

void concater::verbose_files() {
	vrb.title("Files:");

	vrb.bullet("Input sparse files: [" + options["input-sparse-list"].as < string > () + "]");
	vrb.bullet("Output sparse file: [" + options["output-sparse"].as < string > () + "]");

	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void concater::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");
}
