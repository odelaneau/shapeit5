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

#include <phaser/phaser_header.h>

using namespace std;

void phaser::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of thread used");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input", bpo::value< string >(), "Genotypes to be phased")
			("input-region", bpo::value< string >(), "Region to be considered in --input")
			("scaffold", bpo::value< string >(), "Scaffold of haplotypes")
			("scaffold-region", bpo::value< string >(), "Region to be considered in --scaffold")
			("map", bpo::value< string >(), "Genetic map")
			("pedigree", bpo::value< string >(), "Pedigree file (kid father mother")
			("haploids", bpo::value < string >(), "List of haploid samples (e.g. chrX for males). One sample ID per line.");

	bpo::options_description opt_pbwt ("PBWT parameters");
	opt_pbwt.add_options()
			("pbwt-modulo", bpo::value< double >()->default_value(0.1), "Storage frequency of PBWT indexes in cM")
			("pbwt-depth-common", bpo::value< int >()->default_value(2), "Depth of PBWT indexes at common sites to condition on")
			("pbwt-depth-rare", bpo::value< int >()->default_value(2), "Depth of PBWT indexes at rare hets to condition on")
			("pbwt-mac", bpo::value< int >()->default_value(2), "Minimal Minor Allele Count at which PBWT is evaluated")
			("pbwt-mdr", bpo::value < double >()->default_value(0.10), "Maximal Missing Data Rate at which PBWT is evaluated");
	
	bpo::options_description opt_hmm ("HMM parameters");
	opt_hmm.add_options()
			("effective-size", bpo::value<int>()->default_value(15000), "Effective size of the population");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output", bpo::value< string >(), "Phased haplotypes (at common AND rare variants)")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_pbwt).add(opt_hmm).add(opt_output);
}

void phaser::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[SHAPEIT5] Phase_rare (phase rare variants onto a haplotype scaffold)");
	vrb.bullet("Authors       : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 5." + string(PHASE2_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());

	if (args.size() == 0 || options.count("help")) { cout << descriptions << endl; exit(0); }

}

void phaser::check_options() {
	if (!options.count("input")) vrb.error("--input missing");

	if (!options.count("input-region")) vrb.error("--input-region missing");

	if (!options.count("scaffold")) vrb.error("--scaffold missing");

	if (!options.count("scaffold-region")) vrb.error("--scaffold-region missing");

	if (!options.count("output")) vrb.error("--output missing");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");

	if (!options["thread"].defaulted() && !options["seed"].defaulted())
		vrb.warning("Using multi-threading prevents reproducing a run by specifying --seed");

	if (!options["effective-size"].defaulted() && options["effective-size"].as < int > () < 1)
		vrb.error("You must specify a positive effective size");
}

void phaser::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Unphased data : [" + options["input"].as < string > () + "]");
	vrb.bullet("Scaffold data : [" + options["scaffold"].as < string > () + "]");
	if (options.count("map")) vrb.bullet("Genetic Map   : [" + options["map"].as < string > () + "]");
	if (options.count("pedigree")) vrb.bullet("Pedigree file : [" + options["pedigree"].as < string > () + "]");
	vrb.bullet("Output data   : [" + options["output"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void phaser::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");
	vrb.bullet("PBWT    : [depth = " + stb.str(options["pbwt-depth-common"].as < int > ()) + "," + stb.str(options["pbwt-depth-rare"].as < int > ()) + " / modulo = " + stb.str(options["pbwt-modulo"].as < double > ()) + " / mac = " + stb.str(options["pbwt-mac"].as < int > ()) + " / mdr = " + stb.str(options["pbwt-mdr"].as < double > ()) + "]");
	if (options.count("map")) vrb.bullet("HMM     : [Ne = " + stb.str(options["effective-size"].as < int > ()) + " / Recombination rates given by genetic map]");
	else vrb.bullet("HMM     : [Ne = " + stb.str(options["effective-size"].as < int > ()) + " / Constant recombination rate of 1cM per Mb]");
}
