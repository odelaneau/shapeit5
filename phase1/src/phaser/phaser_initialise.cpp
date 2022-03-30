/*******************************************************************************
 * Copyright (C) 2018-2022 Olivier Delaneau, University of Lausanne
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
	if (options.count("scaffold")) readerG.addReferenceFilename(options["scaffold"].as < string > ());
	if (options.count("filter-snp")) readerG.setFilterSNP();
	if (!options["filter-maf"].defaulted()) readerG.setFilterMAF(options["filter-maf"].as < double > ());

	//step2: Read the genotype data
	readerG.scanGenotypes();
	readerG.allocateGenotypes();
	readerG.readGenotypes();

	//step3: Read and initialise genetic map
	vrb.title("Setting up genetic map:");
	if (options.count("map")) {
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < string > ());
		V.setGeneticMap(readerGM);
	} else V.setGeneticMap();
	M.initialise(V, options["hmm-ne"].as < int > (), (readerG.n_main_samples+readerG.n_ref_samples)*2);

	//step4: Initialize haplotype set
	vrb.title("Initializing data structures:");
	G.imputeMonomorphic(V);
	H.updateHaplotypes(G, true);
	H.transposeHaplotypes_H2V(true);
	H.initialize(V,	options["pbwt-modulo"].as < double > (),
					options["pbwt-window"].as < double > (),
					options["pbwt-mdr"].as < double > (),
					options["pbwt-depth"].as < int > (),
					options["pbwt-mac"].as < int > (),
					options["thread"].as < int > ());

	if (!options.count("pbwt-disable-init")) H.solve(&G);

	//step5: Read pedigrees
	if (options.count("pedigree")) {
		pedigree_reader readerP;
		readerP.readPedigreeFile(options["pedigree"].as < string > ());
		G.scaffoldUsingPedigrees(readerP);
	}

	//step6: Initialize genotype structures
	genotype_builder(G, options["thread"].as < int > ()).build();

	//step7: Allocate data structures for computations
	unsigned int max_number_transitions = G.largestNumberOfTransitions();
	unsigned int max_number_missing = G.largestNumberOfMissings();
	threadData = vector < compute_job >(options["thread"].as < int > (), compute_job(V, G, H, max_number_transitions, max_number_missing));
}
