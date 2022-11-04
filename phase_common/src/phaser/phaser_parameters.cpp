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
			("seed", bpo::value < int >()->default_value(15052011), "Seed of the random number generator")
			("thread,T", bpo::value < int >()->default_value(1), "Number of thread used");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,I", bpo::value < string >(), "Genotypes to be phased in VCF/BCF format")
			("reference,H", bpo::value < string >(), "Reference panel of haplotypes in VCF/BCF format")
			("scaffold,S", bpo::value < string >(), "Scaffold of haplotypes in VCF/BCF format")
			("map,M", bpo::value < string >(), "Genetic map")
			("pedigree", bpo::value < string >(), "Pedigree information (kid father mother)")
			("region,R", bpo::value < string >(), "Target region");

	bpo::options_description opt_mcmc ("MCMC parameters");
	opt_mcmc.add_options()
			("mcmc-iterations", bpo::value<string>()->default_value("5b,1p,1b,1p,1b,1p,5m"), "Iteration scheme of the MCMC")
			("mcmc-prune", bpo::value < double >()->default_value(0.999), "Pruning threshold for genotype graphs")
			("mcmc-noinit", "Disable phasing initialization by PBWT sweep");

	bpo::options_description opt_pbwt ("PBWT parameters");
	opt_pbwt.add_options()
			("pbwt-modulo", bpo::value < double >()->default_value(0.1), "Storage frequency of PBWT indexes in cM (i.e. storage every 0.025cM by default)")
			("pbwt-depth", bpo::value < int >()->default_value(4), "Depth of PBWT indexes to condition on")
			("pbwt-mac", bpo::value < int >()->default_value(5), "Minimal Minor Allele Count at which PBWT is evaluated")
			("pbwt-mdr", bpo::value < double >()->default_value(0.1), "Maximal Missing Data Rate at which PBWT is evaluated")
			("pbwt-window", bpo::value < double >()->default_value(4), "Run PBWT selection in windows of this size");
	
	bpo::options_description opt_hmm ("HMM parameters");
	opt_hmm.add_options()
			("hmm-window", bpo::value < double >()->default_value(4), "Minimal size of the phasing window in cM")
			("hmm-ne", bpo::value < int >()->default_value(15000), "Effective size of the population");

	bpo::options_description opt_filter ("FILTER parameters");
	opt_filter.add_options()
			("filter-snp", "Only consider SNPs")
			("filter-maf", bpo::value < double >()->default_value(0.0), "Only consider variants with at lest this MAF, requires AC/AN tags");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Phased haplotypes in VCF/BCF format")
			("output-graph", bpo::value< string >(), "Phased haplotypes in BIN format [Useful to sample multiple likely haplotype configurations per sample]")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_mcmc).add(opt_pbwt).add(opt_hmm).add(opt_filter).add(opt_output);
}

void phaser::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[SHAPEIT5] Phase1 (jointly phase multiple common markers)");
	vrb.bullet("Author        : Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 5." + string(PHASE1_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());
}

void phaser::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("region"))
		vrb.error("You must specify a region or chromosome to phase using --region");

	if ((options.count("output")+options.count("bingraph"))==0)
		vrb.error("You must specify a phased output file with --output");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");

	if (!options["thread"].defaulted() && !options["seed"].defaulted())
		vrb.warning("Using multi-threading prevents reproducing a run by specifying --seed");

	if (!options["hmm-ne"].defaulted() && options["hmm-ne"].as < int > () < 1)
		vrb.error("You must specify a positive effective size");

	if (!options["hmm-window"].defaulted() && (options["hmm-window"].as < double > () < 0.5 || options["hmm-window"].as < double > () > 10))
		vrb.error("You must specify a HMM window size comprised between 0.5 and 10 cM");

	if (!options["pbwt-window"].defaulted() && (options["pbwt-window"].as < double > () < 0.5 || options["pbwt-window"].as < double > () > 10))
		vrb.error("You must specify a PBWT window size comprised between 0.5 and 10 cM");

	parse_iteration_scheme(options["mcmc-iterations"].as < string > ());
}

void phaser::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF     : [" + options["input"].as < string > () + "]");
	if (options.count("reference")) vrb.bullet("Reference VCF : [" + options["reference"].as < string > () + "]");
	if (options.count("scaffold")) vrb.bullet("Scaffold VCF  : [" + options["scaffold"].as < string > () + "]");
	if (options.count("pedigree")) vrb.bullet("Pedigree file : [" + options["pedigree"].as < string > () + "]");
	if (options.count("map")) vrb.bullet("Genetic Map   : [" + options["map"].as < string > () + "]");
	if (options.count("output")) vrb.bullet("Output VCF    : [" + options["output"].as < string > () + "]");
	if (options.count("bingraph")) vrb.bullet("Output BIN    : [" + options["bingraph"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void phaser::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");
	vrb.bullet("MCMC    : " + get_iteration_scheme());

	pbwt_auto = options["pbwt-modulo"].defaulted() && options["pbwt-depth"].defaulted();
	if (!pbwt_auto)
		vrb.bullet("PBWT    : [window = " + stb.str(options["pbwt-window"].as < double > ()) + "cM / depth = " + stb.str(options["pbwt-depth"].as < int > ()) + " / modulo = " + stb.str(options["pbwt-modulo"].as < double > ()) + " / mac = " + stb.str(options["pbwt-mac"].as < int > ()) + " / missing = " + stb.str(options["pbwt-mdr"].as < double > ()) + "]");
	else
		vrb.bullet("PBWT    : [window = " + stb.str(options["pbwt-window"].as < double > ()) + "cM / depth = auto / modulo = auto / mac = " + stb.str(options["pbwt-mac"].as < int > ()) + " / missing = " + stb.str(options["pbwt-mdr"].as < double > ()) + "]");

	if (options.count("map"))  vrb.bullet("HMM     : [window = " + stb.str(options["hmm-window"].as < double > ()) + "cM / Ne = " + stb.str(options["hmm-ne"].as < int > ()) + " / Recombination rates given by genetic map]");
	else vrb.bullet("HMM     : [window = " + stb.str(options["hmm-window"].as < double > ()) + "cM / Ne = " + stb.str(options["hmm-ne"].as < int > ()) + " / Constant recombination rate of 1cM per Mb]");
	if (options.count("filter-snp") || (!options["filter-maf"].defaulted()))
		vrb.bullet("FILTERS : [snp only = " + stb.str(options.count("filter-snp")) + " / MAF = " + stb.str(options["filter-maf"].as < double > ()) + "]");
}
