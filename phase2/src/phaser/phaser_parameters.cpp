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

#include "../../versions/versions.h"

#include <phaser/phaser_header.h>

void phaser::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of thread used");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input", bpo::value< string >(), "Genotypes to be phased in VCF/BCF format")
			("input-region", bpo::value< string >(), "Region to be considered in --input")
			("input-maf", bpo::value< double >()->default_value(0.001), "Threshold for sparse genotype representation")
			("scaffold", bpo::value< string >(), "Scaffold of haplotypes in VCF/BCF format")
			("scaffold-region", bpo::value< string >(), "Region to be considered in --scaffold")
			("map", bpo::value< string >(), "Genetic map");
/*
	bpo::options_description opt_sequence ("Sequencing data files");
	opt_input.add_options()
			("bam-list", bpo::value< string >(), "List of BAM files to be used for read based phasing")
			("bam-fasta", bpo::value< string >(), "Fasta file for CRAM decompression")
			("bam-mapq", bpo::value< int >()->default_value(10), "Minimal mapping quality to consider a read")
			("bam-baseq", bpo::value< int >()->default_value(20), "Minimal calling quality to consider a base in a read");
*/
	bpo::options_description opt_mcmc ("MCMC parameters");
	opt_mcmc.add_options()
			("mcmc-iterations", bpo::value< int >()->default_value(10), "Number of MCMC iterations")
			("mcmc-burnin", bpo::value< int >()->default_value(5), "Number of MCMC burn-in iterations");

	bpo::options_description opt_pbwt ("PBWT parameters");
	opt_pbwt.add_options()
			("pbwt-modulo", bpo::value< double >()->default_value(0.1), "Storage frequency of PBWT indexes in cM")
			("pbwt-depth-common", bpo::value< int >()->default_value(4), "Depth of PBWT indexes at common sites to condition on")
			("pbwt-depth-rare", bpo::value< int >()->default_value(4), "Depth of PBWT indexes at rare hets to condition on")
			("pbwt-mac", bpo::value< int >()->default_value(2), "Minimal Minor Allele Count at which PBWT is evaluated")
			("pbwt-mdr", bpo::value < double >()->default_value(0.10), "Maximal Missing Data Rate at which PBWT is evaluated");
	
	bpo::options_description opt_hmm ("HMM parameters");
	opt_hmm.add_options()
			("effective-size", bpo::value<int>()->default_value(15000), "Effective size of the population")
			("compress-threshold", bpo::value<double>()->default_value(0.01), "Ignore states with probability below this threshold");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Phased haplotypes in VCF/BCF format")
			("log", bpo::value< string >(), "Log file");

	//descriptions.add(opt_base).add(opt_input).add(opt_sequence).add(opt_mcmc).add(opt_pbwt).add(opt_hmm).add(opt_output);
	descriptions.add(opt_base).add(opt_input).add(opt_mcmc).add(opt_pbwt).add(opt_hmm).add(opt_output);
}

void phaser::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[SHAPEIT5] Phase rare variants and indels onto a dense haplotype scaffold");
	vrb.bullet("Authors       : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@gmail.com");
	vrb.bullet("Version       : " + string(PHASE2_VERSION));
	vrb.bullet("Run date      : " + tac.date());
}

void phaser::check_options() {
	if (!options.count("input"))
		vrb.error("--input missing");

	if (!options.count("input-region"))
		vrb.error("--input-region missing");

	if (!options.count("scaffold"))
		vrb.error("--scaffold missing");

	if (!options.count("scaffold-region"))
		vrb.error("--scaffold-region missing");

	if (!options.count("output"))
		vrb.error("You must specify a phased output file with --output");

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
	vrb.bullet("Input VCF     : [" + options["input"].as < string > () + "]");
	vrb.bullet("Scaffold VCF  : [" + options["scaffold"].as < string > () + "]");
	if (options.count("map")) vrb.bullet("Genetic Map   : [" + options["map"].as < string > () + "]");
	//if (options.count("bam-list")) vrb.bullet("BAM list      : [" + options["bam-list"].as < string > () + "]");
	if (options.count("output")) vrb.bullet("Output VCF    : [" + options["output"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void phaser::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");
	vrb.bullet("MCMC    : #iterations=" + stb.str(options["mcmc-iterations"].as < int > ()) + " / #burn-in=" + stb.str(options["mcmc-burnin"].as < int > ()));
	vrb.bullet("PBWT    : [depth = " + stb.str(options["pbwt-depth-common"].as < int > ()) + "," + stb.str(options["pbwt-depth-rare"].as < int > ()) + " / modulo = " + stb.str(options["pbwt-modulo"].as < double > ()) + " / mac = " + stb.str(options["pbwt-mac"].as < int > ()) + " / mdr = " + stb.str(options["pbwt-mdr"].as < double > ()) + "]");
	if (options.count("map")) vrb.bullet("HMM     : [Ne = " + stb.str(options["effective-size"].as < int > ()) + " / Recombination rates given by genetic map]");
	else vrb.bullet("HMM     : [Ne = " + stb.str(options["effective-size"].as < int > ()) + " / Constant recombination rate of 1cM per Mb]");
	//if (options.count("bam-list")) vrb.bullet("BAM     : [MappingQ = " + stb.str(options["bam-mapq"].as < int > ()) + " / BaseQ = " + stb.str(options["bam-baseq"].as < int > ()) + "]");
}
