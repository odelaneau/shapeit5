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

#include <phaser/phaser_header.h>

#include <io/genotype_reader/genotype_reader_header.h>
#include <io/haplotype_writer.h>
#include <io/gmap_reader.h>
#include <io/pedigree_reader.h>
#include <modules/genotype_builder.h>

using namespace std;

void phaser::read_files_and_initialise() {
	//step0: Initialize seed and multi-threading
	rng.setSeed(options["seed"].as < int > ());
	if (options["thread"].as < int > () > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//step1: Set up the genotype reader
	vrb.title("Reading genotype data:");
	genotype_reader readerG(H, G, V);
	readerG.setThreads(options["thread"].as < int > ());
	readerG.setRegion(options["region"].as < string > ());
	readerG.setMainFilename(options["input"].as < string > ());
	if (options.count("reference")) readerG.addReferenceFilename(options["reference"].as < string > ());
	if (options.count("scaffold")) readerG.addScaffoldFilename(options["scaffold"].as < string > ());
	if (options.count("filter-snp")) readerG.setFilterSNP();
	if (!options["filter-maf"].defaulted()) readerG.setFilterMAF(options["filter-maf"].as < double > ());

	//step2: Read the genotype data
	readerG.scanGenotypes();
	readerG.allocateGenotypes();
	readerG.readGenotypes();

	//step3: Read pedigrees
	if (options.count("pedigree")) {
		pedigree_reader readerP;
		readerP.readPedigreeFile(options["pedigree"].as < string > ());
		G.scaffoldUsingPedigrees(readerP);
	}

	//step4: Read and initialise genetic map
	vrb.title("Setting up genetic map:");
	if (options.count("map")) {
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < string > ());
		V.setGeneticMap(readerGM);
	} else V.setGeneticMap();
	M.initialise(V, options["hmm-ne"].as < int > (), (readerG.n_main_samples+readerG.n_ref_samples)*2);

	//step5: Initialize haplotype set
	vrb.title("Initializing data structures:");
	G.imputeMonomorphic(V);
	H.updateHaplotypes(G, true);
	H.transposeHaplotypes_H2V(true);

	//step6: Initialize PBWT for selecting states
	if (pbwt_auto) {
		unsigned int cumulative_sample_size = readerG.n_main_samples + readerG.n_ref_samples;
		pbwt_depth = max(min((int)round(10-log10(cumulative_sample_size)), 8), 2);
		pbwt_modulo = max(min((log(cumulative_sample_size) - log(50) + 1) * 0.01, 0.15), 0.005);
		vrb.bullet("PBWT parameters auto setting : [modulo = " + stb.str(pbwt_modulo, 3) + " / depth = " + stb.str(pbwt_depth, 3) + "]");
	} else {
		pbwt_depth = options["pbwt-depth"].as < int > ();
		pbwt_modulo = options["pbwt-modulo"].as < double > ();
	}

	H.initialize(V,	pbwt_modulo,
					options["pbwt-window"].as < double > (),
					options["pbwt-mdr"].as < double > (),
					pbwt_depth,
					options["pbwt-mac"].as < int > (),
					options["thread"].as < int > ());

	if (!options.count("pbwt-disable-init")) H.solve(&G);

	//step6: Initialize genotype structures
	genotype_builder(G, options["thread"].as < int > ()).build();

	//step7: Allocate data structures for computations
	unsigned int max_number_transitions = G.largestNumberOfTransitions();
	unsigned int max_number_missing = G.largestNumberOfMissings();
	threadData = vector < compute_job >(options["thread"].as < int > (), compute_job(V, G, H, max_number_transitions, max_number_missing));
}
