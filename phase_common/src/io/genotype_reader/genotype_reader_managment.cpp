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

#include <io/genotype_reader/genotype_reader_header.h>

using namespace std;

genotype_reader::genotype_reader(haplotype_set & _H, genotype_set & _G, variant_map & _V) : H(_H), G(_G), V(_V) {
	nthreads = 1;
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
	region = "";
	n_genotypes = vector < unsigned long > (5, 0);
	filter_min_maf = -1.0f;
	filter_snp_only = false;
	filenames = vector < string > (3, "");
	panels = vector < char > (3, 0) ;
}

genotype_reader::~genotype_reader() {
	nthreads = 1;
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
	region = "";
	n_genotypes = vector < unsigned long > (5, 0);
	filter_min_maf = -1.0f;
	filter_snp_only = false;
	filenames = vector < string > (3, "");
	panels = vector < char > (3, 0) ;
}

void genotype_reader::allocateGenotypes() {
	assert(n_variants != 0 && (n_main_samples+n_ref_samples) != 0);
	G.allocate(n_main_samples, n_variants);
	H.allocate(n_main_samples, n_ref_samples, n_variants);
}

void genotype_reader::setFilterMAF(float _filter_min_maf) { filter_min_maf = _filter_min_maf; }

void genotype_reader::setFilterSNP() { filter_snp_only = true; }

void genotype_reader::setMainFilename(string file) { filenames[0] = file; panels[0] = 1; }

void genotype_reader::addReferenceFilename(string file) { filenames[1] = file; panels[1] = 1; }

void genotype_reader::addScaffoldFilename(string file) { filenames[2] = file; panels[2] = 1; }

void genotype_reader::setThreads(int _threads) { threads = _threads; }

void genotype_reader::setRegion(string _region) { region = _region; }

