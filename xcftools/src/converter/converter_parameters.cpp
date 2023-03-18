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
			("thread,T", bpo::value<int>()->default_value(1), "Number of thread used for VCF/BCF (de-)compression");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,i", bpo::value< string >(), "Input genotype data in plain VCF/BCF format")
			("region,r", bpo::value< string >(), "Region to be considered in --input")
			("maf,m", bpo::value< float >()->default_value(0.001), "Threshold to distinguish rare variants from common ones");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,o", bpo::value< string >()->default_value("-"), "Output file [- for stdout]")
			("format,O", bpo::value< string >()->default_value("bcf"), "Output file format")
			("log", bpo::value< string >(), "Output log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void converter::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	string output = options["output"].as < string > ();
	if (output == "-") vrb.set_silent();

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[XCFtools] Convert from/to XCF files");
	vrb.bullet("Authors       : Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 5." + string(XCFTLS_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());
}

void converter::check_options() {
	if (!options.count("input")) vrb.error("--input needs to be specified");
	if (!options.count("region")) vrb.error("--region needs to be specified");
	if (!options.count("format")) vrb.error("--format needs to be specified");

	string format = options["format"].as < string > ();
	string output = options["output"].as < string > ();
	if (!isBCF(format) && output == "-") vrb.error("Only BCF format [bcf] is supported on stdout");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");
}

void converter::verbose_files() {
	vrb.title("Files:");
	string format = options["format"].as < string > ();
	string output = options["output"].as < string > ();
	if (isXCF(format)) {
		vrb.bullet("Input BCF     : [" + options["input"].as < string > () + "]");
		vrb.bullet("Output XCF    : [" + output + "]");
	} else if (isBCF(format)) {
		vrb.bullet("Input XCF     : [" + options["input"].as < string > () + "]");
		if (output == "-") vrb.bullet("Output BCF   : [STDOUT] / uncompressed");
		else vrb.bullet("Output BCF    : [" + options["output"].as < string > () + "]");
	} else vrb.error("Output format [" + format + "] unrecognized");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void converter::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");

	string format = options["format"].as < string > ();
	if (format[0] == 's') vrb.bullet("MAF     : " + stb.str(options["maf"].as < float > ()));
}
