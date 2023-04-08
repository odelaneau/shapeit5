/*******************************************************************************
 * Copyright (C) 2020 Olivier Delaneau, University of Lausanne
 * Copyright (C) 2020 Simone Rubinacci, University of Lausanne
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

#include <family/family_header.h>


void family::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input", bpo::value < std::string >(), "Genotypes to be phased in VCF/BCF format")
			("input-haploid", "Input data is haploid (1 allele per sample)")
			("array", bpo::value < std::string >(), "Positions in a SNP array")
			("region", bpo::value < std::string >(), "Whole genomic region to be phased (including left/right buffers)");

	bpo::options_description opt_algo ("Parameters");
	opt_algo.add_options()
			("missing", bpo::value< double >()->default_value(0.0), "Fraction of missing genotypes")
			("twins", bpo::value< int >()->default_value(0), "Number of twins to include")
			("obs-fam", bpo::value< int >()->default_value(0), "Number of full families to include")
			("obs-kid", bpo::value< int >()->default_value(0), "Number of kids per full family included")
			("hid-fam", bpo::value< int >()->default_value(0), "Number of hidden families to include")
			("hid-kid", bpo::value< int >()->default_value(0), "Number of kids per hidden family included");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("validation", bpo::value< std::string >(), "VCF/BCF file reading for validation")
			("estimation", bpo::value< std::string >(), "VCF/BCF file reading for phasing")
			("log", bpo::value< std::string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_algo).add(opt_output);
}

void family::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < std::string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < std::string > () +"]");

	vrb.title("[SIMULATE] Simulate family haplotypes data");
	vrb.bullet("Author        : Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@unil.ch");
	vrb.bullet("Run date      : " + tac.date());

	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }
}

void family::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("region"))
		vrb.error("You must specify a region to phase using --region");

	if (!options.count("estimation"))
		vrb.error("You must specify a region to output using --estimation");

	if (!options.count("validation"))
		vrb.error("You must specify a region to output using --validation");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");
}

void family::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF     : [" + options["input"].as < std::string > () + "]");
	vrb.bullet("Validation VCF: [" + options["validation"].as < std::string > () + "]");
	vrb.bullet("Estimation VCF: [" + options["estimation"].as < std::string > () + "]");
	vrb.bullet("Input region  : [" + options["region"].as < std::string > () + "]");
	if (options.count("array")) vrb.bullet("Array FILE    : [" + options["array"].as < std::string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < std::string > () + "]");
}

void family::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed       : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Haploid inp: " + stb.str(options.count("input-haploid")));
	vrb.bullet("Missing    : " + stb.str(options["missing"].as < double > ()));
	vrb.bullet("Twins      : " + stb.str(options["twins"].as < int > ()));
	vrb.bullet("Full fam.  : " + stb.str(options["obs-fam"].as < int > ()));
	vrb.bullet("Full kid.  : " + stb.str(options["obs-kid"].as < int > ()));
	vrb.bullet("Hidden fam.: " + stb.str(options["hid-fam"].as < int > ()));
	vrb.bullet("Hidden kid.: " + stb.str(options["hid-kid"].as < int > ()));
}
