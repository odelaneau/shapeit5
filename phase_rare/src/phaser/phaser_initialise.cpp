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

using namespace std;

void phaser::read_files_and_initialise() {
	//step0: Initialize seed
	rng.setSeed(options["seed"].as < int > ());
	nthreads = options["thread"].as < int > ();
	if (nthreads > 1) {
		i_jobs = 0;
		id_workers = vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

    //step1: Parsing region string
	buildCoordinates();

	//step2: Set up the genotype reader
	vrb.title("Reading genotype data:");
	genotype_reader readerG(H, G, V);
	readerG.setThreads(options["thread"].as < int > ());
	readerG.setRegions(options["scaffold-region"].as < string > (), input_start, input_stop);
	readerG.setMAF(options["input-maf"].as < double > ());

	//step3: Read the genotype data
	if (options.count("input-plain")) {
		readerG.setFilenames(options["input-plain"].as < string > (), options["input-plain"].as < string > (), options["scaffold"].as < string > ());
		readerG.scanGenotypesPlain();
		readerG.allocateGenotypes();
		readerG.readGenotypesPlain();
	} else {
		readerG.setFilenames(options["input-sparse"].as < string > () + ".bcf", options["input-sparse"].as < string > () + ".bin", options["scaffold"].as < string > ());
		readerG.scanGenotypesSparse();
		readerG.allocateGenotypes();
		readerG.readGenotypesSparse();
	}


	G.imputeMonomorphic();
	G.fillup_by_transpose_V2I();

	//step4: Read pedigrees
	if (options.count("pedigree")) {
		pedigree_reader readerP;
		readerP.readPedigreeFile(options["pedigree"].as < string > ());
		G.phaseUsingPedigrees(readerP);
	}

	//step5: Read and initialise genetic map
	vrb.title("Setting up genetic map:");
	if (options.count("map")) {
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < string > ());
		V.setGeneticMap(readerGM);
	} else V.setGeneticMap();
	M.initialise(V, options["effective-size"].as < int > (), readerG.n_samples*2);
	/*
	double theta = 1.0f / (log(readerG.n_samples*2) + 0.5);
	rare_genotype::ed = theta / (2*( readerG.n_samples*2 + theta ));
	rare_genotype::ee = 1.0f - rare_genotype::ed;
	vrb.bullet("Emission probability = " + stb.str(rare_genotype::ed));
	*/
	H.transposeHaplotypes_V2H();
}
