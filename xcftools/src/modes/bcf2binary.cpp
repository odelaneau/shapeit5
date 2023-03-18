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

#include <modes/bcf2binary.h>
#include <utils/xcf.h>

#include <containers/bitvector.h>
#include <objects/rare_genotype.h>

using namespace std;

bcf2binary::bcf2binary(string _region, float _minmaf, int _nthreads, int _mode) {
	mode = _mode;
	nthreads = _nthreads;
	region = _region;
	minmaf = _minmaf;
	vector < string > tokens;
	stb.split(region, tokens, ":");
	if (tokens.size() == 1 || tokens.size() == 2) contig = tokens[0];
	else vrb.error("Could not parse --region");
}

bcf2binary::~bcf2binary() {
}

void bcf2binary::convert(string finput, string foutput) {
	tac.clock();

	switch (mode) {
	case CONV_BCF_BG: vrb.title("Converting from BCF to XCF [Binary/Genotype]"); break;
	case CONV_BCF_BH: vrb.title("Converting from BCF to XCF [Binary/Haplotype]"); break;
	case CONV_BCF_SG: vrb.title("Converting from BCF to XCF [Sparse/Genotype]"); break;
	case CONV_BCF_SH: vrb.title("Converting from BCF to XCF [Sparse/Haplotype]"); break;
	}
	vrb.bullet("Region        : " + region);
	vrb.bullet("Contig        : " + contig);
	if (mode == CONV_BCF_SG || mode == CONV_BCF_SH) vrb.bullet("Min MAF       : " + stb.str(minmaf));

	//Opening XCF reader for input
	xcf_reader XR(region, nthreads);
	int32_t idx_file = XR.addFile(finput);

	//Check file type
	int32_t type = XR.typeFile(idx_file);
	if (type != FILE_BCF) vrb.error("[" + finput + "] is not a BCF file");

	//Get sample IDs
	vector < string > samples;
	int32_t nsamples = XR.getSamples(idx_file, samples);
	vrb.bullet("#samples = " + stb.str(nsamples));

	//Opening XCF writer for output [false means NO records in BCF body but in external BIN file]
	xcf_writer XW(foutput, false, nthreads);

	//Write header
	XW.writeHeader(samples, contig, string("XCFtools ") + string(XCFTLS_VERSION));

	//Allocate input/output buffer
	int32_t * input_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));
	int32_t * output_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));

	//Allocate output bitvector for common variants
	bitvector binary_buffer = bitvector (2 * nsamples);

	//Proceed with conversion
	uint32_t n_lines_rare = 0, n_lines_comm = 0;
	while (XR.nextRecord()) {

		//Copy over variant information
		XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());

		//Is that a rare variant?
		float af =  XR.getAF();
		float maf = min(af, 1.0f-af);
		bool minor = (af < 0.5f);
		bool rare = (maf < minmaf);

		//Get record
		XR.readRecord(0, reinterpret_cast< char** > (&input_buffer));

		//Convert
		uint32_t n_sparse = 0;
		for (uint32_t i = 0 ; i < nsamples ; i++) {
			bool a0 = (bcf_gt_allele(input_buffer[2*i+0])==1);
			bool a1 = (bcf_gt_allele(input_buffer[2*i+1])==1);
			bool mi = (input_buffer[i+0] == bcf_gt_missing || input_buffer[i+1] == bcf_gt_missing);

			if (mi && (mode == CONV_BCF_SH || mode == CONV_BCF_BH)) vrb.error("Missing data in phased data is not permitted!");

			//BCF => SPARSE GENOTYPE
			if (mode == CONV_BCF_SG) {
				if (rare) {
					if (a0 == minor || a1 == minor || mi)
						output_buffer[n_sparse++] = rare_genotype(i, (a0!=a1), mi, a0, a1, 0).get();
				} else {
					if (mi) { binary_buffer.set(2*i+0, true); binary_buffer.set(2*i+1, false); }		//Missing as 10
					else if (a0 == a1) { binary_buffer.set(2*i+0, a0); binary_buffer.set(2*i+1, a1); }
					else { binary_buffer.set(2*i+0, false); binary_buffer.set(2*i+1, true); }			//Hets as 01
				}
			}

			//BCF => SPARSE HAPLOTYPE
			if (mode == CONV_BCF_SH) {
				if (rare) {
					if (a0 == minor) output_buffer[n_sparse++] = 2*i+0;
					if (a1 == minor) output_buffer[n_sparse++] = 2*i+1;
				} else {
					binary_buffer.set(2*i+0, a0);
					binary_buffer.set(2*i+1, a1);
				}
			}

			//BCF => BINARY GENOTYPE
			if (mode == CONV_BCF_BG) {
				if (mi) { binary_buffer.set(2*i+0, true); binary_buffer.set(2*i+1, false); }		//Missing as 10
				else if (a0 == a1) { binary_buffer.set(2*i+0, a0); binary_buffer.set(2*i+1, a1); }
				else { binary_buffer.set(2*i+0, false); binary_buffer.set(2*i+1, true); }			//Hets as 01
			}

			//BCF => BINARY HAPLOTYPE
			if (mode == CONV_BCF_BH) {
				binary_buffer.set(2*i+0, a0);
				binary_buffer.set(2*i+1, a1);
			}
		}

		//Write record
		if (mode == CONV_BCF_SG && rare)
			XW.writeRecord(RECORD_SPARSE_GENOTYPE, reinterpret_cast<char*>(output_buffer), n_sparse * sizeof(int32_t));
		else if (mode == CONV_BCF_SH && rare)
			XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(output_buffer), n_sparse * sizeof(int32_t));
		else if (mode == CONV_BCF_SG || mode == CONV_BCF_BG)
			XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_buffer.bytes, binary_buffer.n_bytes);
		else
			XW.writeRecord(RECORD_BINARY_HAPLOTYPE, binary_buffer.bytes, binary_buffer.n_bytes);

		//Line counting
		n_lines_comm += !rare || mode == CONV_BCF_BG || mode == CONV_BCF_BH;
		n_lines_rare += rare && (mode == CONV_BCF_SG || mode == CONV_BCF_SH);

		//Verbose
		if ((n_lines_comm+n_lines_rare) % 10000 == 0) {
			if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of BCF records processed: N=" + stb.str(n_lines_comm));
			else vrb.bullet("Number of BCF records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));
		}
	}

	if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of BCF records processed: N=" + stb.str(n_lines_comm));
	else vrb.bullet("Number of BCF records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));

	//Free
	free(input_buffer);
	free(output_buffer);

	//Close files
	XR.close();
	XW.close();
}

