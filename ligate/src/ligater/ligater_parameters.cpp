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

#include "../../versions/versions.h"
#include <ligater/ligater_header.h>

using namespace std;

void ligater::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input", bpo::value < string >(), "Text file containing all VCF/BCF to ligate")
			("pedigree", bpo::value< string >(), "Pedigree file (kid father mother");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Output ligated file in VCF/BCF format")
			("index", "Whether to index the ligated output (csi format)")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void ligater::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	vrb.title("[SHAPEIT5] Ligate (ligate multiple output files into chromosome-wide files)");
	vrb.bullet("Authors       : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 5." + string(LIGATE_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");
}

void ligater::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify the list of VCF/BCF to ligate using --input");

	if (!options.count("output"))
		vrb.error("You must specify an output VCF/BCF file with --output");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options["thread"].as < int > () < 1)
		vrb.error("Number of threads is a strictly positive number.");
}

void ligater::verbose_files() {
	std::array<std::string,2> no_yes = {"NO","YES"};

	vrb.title("Files:");
	vrb.bullet("Input LIST     : [" + options["input"].as < string > () + "]");
	vrb.bullet("Output VCF     : [" + options["output"].as < string > () + "]");
	vrb.bullet("Index output   : [" + no_yes[options.count("index")] + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void ligater::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed           : [" + stb.str(options["seed"].as < int > ()) + "]");
	vrb.bullet("#Threads       : [" + stb.str(options["thread"].as < int > ()) + "]");
}
