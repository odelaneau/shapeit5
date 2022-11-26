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


#include <switcher/switcher_header.h>

using namespace std;

void switcher::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("thread,T", bpo::value<int>()->default_value(1), "Number of thread used");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("validation,V", bpo::value< string >(), "Validation dataset in VCF/BCF format")
			("estimation,E", bpo::value< string >(), "Phased dataset in VCF/BCF format")
			("frequency,F", bpo::value< string >(), "Variant frequency in VCF/BCF format")
			("pedigree,P", bpo::value< string >(), "Pedigree file in PED format")
			("region,R", bpo::value< string >(), "Target region")
			("nbins", bpo::value<int>()->default_value(20), "Number of bins used for calibration")
			("min-pp", bpo::value<double>()->default_value(0.0f), "Minimal PP value for entering computations")
			("singleton", "Trick to phase singleton")
			("dupid", "Duplicate ID for UKB matching IDs");
	
	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Prefix for all report files")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void switcher::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[SHAPEIT5] Switch (validation of haplotypes using trios/duos or true sets)");
	vrb.bullet("Author        : Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 5." + string(SWITCH_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());

}

void switcher::check_options() {
	if (!options.count("validation")) vrb.error("You must specify --validation");
	if (!options.count("estimation")) vrb.error("You must specify --estimation");
	if (!options.count("region")) vrb.error("You must specify a region or chromosome to process using --region");
	if (!options.count("output")) vrb.error("You must specify a prefix for output files with --output");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");
}

void switcher::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Validation VCF: [" + options["validation"].as < string > () + "]");
	vrb.bullet("Phased VCF    : [" + options["estimation"].as < string > () + "]");
	if (options.count("frequency")) vrb.bullet("Frequency VCF : [" + options["frequency"].as < string > () + "]");
	vrb.bullet("Output prefix : [" + options["output"].as < string > () + "]");
	if (options.count("pedigree")) vrb.bullet("Pedigree file : [" + options["pedigree"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void switcher::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("#threads : " + stb.str(options["thread"].as < int > ()));
	vrb.bullet("#bins    : " + stb.str(options["nbins"].as < int > ()));
	vrb.bullet("MinPP    : " + stb.str(options["min-pp"].as < double > ()));
}
