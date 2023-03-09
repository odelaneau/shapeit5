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

#include <converter/converter_header.h>

using namespace std;

void converter::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of thread used for VCF/BCF (de-)compression");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input-plain", bpo::value< string >(), "Input genotype data in plain VCF/BCF format")
			("input-sparse", bpo::value< vector < string > >()->multitoken(), "Input genotype data in sparse VCF/BCF format")
			("region", bpo::value< string >(), "Region to be considered in --input")
			("maf", bpo::value< double >()->default_value(0.001), "Threshold for sparse genotype representation");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output-plain", bpo::value< string >(), "Output genotype data in plain VCF/BCF format")
			("output-sparse", bpo::value< vector < string > >()->multitoken(), "Output genotype data in sparse VCF/BCF format")
			("force-unphased", "Force out sparse to be unphased")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void converter::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[SHAPEIT5] Convert from/to sparse VCF/BCF");
	vrb.bullet("Authors       : Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 5." + string(SCFTLS_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());
}

void converter::check_options() {
	if (options.count("input-sparse")) {
		if (!options.count("output-plain")) vrb.error("Use --input-sparse together with --output-plain");
		vector < string > files = options["input-sparse"].as < vector < string > > ();
		if (files.size() != 2) vrb.error("--input-sparse needs 2 arguments");
	} else if (options.count("input-plain")) {
		if (!options.count("output-sparse")) vrb.error("Use --input-plain together with --output-sparse");
		vector < string > files = options["output-sparse"].as < vector < string > > ();
		if (files.size() != 2) vrb.error("--output-sparse needs 2 arguments");
	} else vrb.error("Use [--input-sparse / --output-plain] OR [--input-plain / --output-sparse]");

	if (!options.count("region"))
		vrb.error("--region missing");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");
}

void converter::verbose_files() {
	vrb.title("Files:");

	int mode1 = options.count("input-plain") + options.count("output-sparse");
	int mode2 = options.count("input-sparse") + options.count("output-plain");

	if (mode1==2 && mode2==0) {
		vector < string > files = options["output-sparse"].as < vector < string > > ();
		vrb.bullet("Input plain VCF/BCF : [" + options["input-plain"].as < string > () + "]");
		vrb.bullet("Output sparse files: [" + files[0] + " / " + files[1] + "]");
	} else if (mode1==0 && mode2==2) {
		vector < string > files = options["input-sparse"].as < vector < string > > ();
		vrb.bullet("Input sparse prefix : [" + files[0] + " / " + files[1] + "]");
		vrb.bullet("Output plain VCF/BCF: [" + options["output-plain"].as < string > () + "]");
	} else {
		vrb.error("Unrecognized conversion mode. Use --input-plain with --output-sparse OR --input-sparse with --output-plain.");
	}

	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void converter::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");
	int mode1 = options.count("input-plain") + options.count("output-sparse");
	int mode2 = options.count("input-sparse") + options.count("output-plain");
	if (mode1 == 2) vrb.bullet("Mode    : Plain VCF to sparse VCF");
	if (mode2 == 2) vrb.bullet("Mode    : Sparse VCF to plain VCF");
}
