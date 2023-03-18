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

#include <utils/xcf.h>
#include <io/genotype_reader/genotype_reader_header.h>

using std::vector;
using std::min;

void genotype_reader::scanGenotypes() {
	tac.clock();
	vrb.wait("  * BCF scanning");

	//Initialize XCF reader
	xcf_reader XR (scaffold_region, nthreads);

	//Open files
	uint32_t idx_unphased = XR.addFile(funphased);
	uint32_t idx_scaffold = XR.addFile(fphased);

	//Sample processing
	std::vector < std::string > samples_unphased, samples_scaffold;
	n_samples = XR.getSamples(idx_unphased, samples_unphased);
	XR.getSamples(idx_scaffold, samples_scaffold);
	if (samples_scaffold != samples_unphased) vrb.error("Different sets of samples across files");
	else G.names = samples_unphased;

	//Parsing BCFs
	uint32_t n_multi = 0, n_both = 0, n_buffer = 0;
	while (XR.nextRecord()) {

		//In Scaffold
		if (XR.hasRecord(idx_scaffold)) {
			V.push(new variant (XR.chr, XR.pos, XR.rsid, XR.ref, XR.alt, (XR.getAF(idx_scaffold) < 0.5f), VARTYPE_SCAF));
			n_scaffold_variants++;
			n_total_variants ++;
			n_both += XR.hasRecord(idx_unphased);
		}

		//In Unphased
		else if (XR.hasRecord(idx_unphased)) {
			if (XR.pos >= input_start && XR.pos <= input_stop) {
				V.push(new variant (XR.chr, XR.pos, XR.rsid, XR.ref, XR.alt, (XR.getAF(idx_unphased) < 0.5f), VARTYPE_RARE));
				n_rare_variants ++; n_total_variants ++;
			} else n_buffer++;
		}

		//Not in both files: multi-allelic variant
		else n_multi ++;
	}

	//Close
	XR.close();

	//Verbose
	if (n_rare_variants == 0) vrb.error("No variants to be phased!");
	vrb.bullet("BCF scanning (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("Variant sites: [#scaffold=" + stb.str(n_rare_variants) + " | #rare=" + stb.str(n_scaffold_variants) + " | #all=" + stb.str(n_total_variants) + "]");
	if (n_buffer) vrb.bullet2(stb.str(n_buffer) + " rare variants in buffer regions excluded");
	if (n_both) vrb.bullet2(stb.str(n_both) + " variants found in both files [priority to scaffold]");
	if (n_multi) vrb.bullet2(stb.str(n_multi) + " multi-allelic variants excluded");
}
