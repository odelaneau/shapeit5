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

genotype_reader::genotype_reader(haplotype_set & _H, genotype_set & _G, variant_map & _V) : H(_H), V(_V), G(_G) {
	nthreads = 1;
	n_scaffold_variants = 0;
	n_rare_variants = 0;
	n_samples = 0;
	funphased = "";
	fphased = "";
	foutput = "";
	scaffold_region = "";
	input_start = 0;
	input_start = 1000000000;
	n_scaffold_genotypes = vector < unsigned long > (4, 0);
	n_rare_genotypes = vector < unsigned long > (4, 0);
}

genotype_reader::~genotype_reader() {
	nthreads = 1;
	n_scaffold_variants = 0;
	n_rare_variants = 0;
	n_samples = 0;
	funphased = "";
	fphased = "";
	scaffold_region = "";
	input_start = 0;
	input_start = 1000000000;
	n_scaffold_genotypes = vector < unsigned long > (4, 0);
	n_rare_genotypes = vector < unsigned long > (4, 0);
}

void genotype_reader::allocateGenotypes() {
	H.allocate(n_samples, n_scaffold_variants);
	G.allocate(V, n_samples, n_scaffold_variants, n_rare_variants);
}

void genotype_reader::setFilenames (string _funphased, string _fbinary, string _fphased, string _foutput) {
	fphased = _fphased;
	fbinary = _fbinary;
	funphased = _funphased;
	foutput = _foutput;
}

void genotype_reader::setThreads(int _nthreads) { nthreads = _nthreads; }

void genotype_reader::setRegions(string _scaffold_region, int _input_start, int _input_stop) {
	input_start = _input_start;
	input_stop = _input_stop;
	scaffold_region = _scaffold_region;
}

void genotype_reader::setMAF(float _minmaf) { minmaf = _minmaf; }
