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

#include <io/haplotype_writer.h>
#include <utils/xcf.h>
#include <containers/bitvector.h>

using namespace std;


haplotype_writer::haplotype_writer(haplotype_set & _H, genotype_set & _G, variant_map & _V, uint32_t _nthreads): H(_H), G(_G), V(_V) {
	nthreads = _nthreads;
}

haplotype_writer::~haplotype_writer() {
}

void haplotype_writer::writeHaplotypes(string fname, string fformat) {
	tac.clock();

	//Open XCF writer
	bool hts_genotypes = (fformat == "bcf");
	xcf_writer XW(fname, hts_genotypes, nthreads);

	//Write header
	vector < string > snames;
	for (int32_t i = 0 ; i < G.n_ind ; i ++) snames.push_back(G.vecG[i]->name.c_str());
	XW.writeHeader(snames, V.vec_pos[0]->chr, string("SHAPEIT5 phase_common ") + string(PHASE1_VERSION));

	//Allocate buffers
	int32_t * output_buffer = (int32_t *) malloc(G.n_ind * 2 * sizeof(int32_t));
	bitvector output_bitvector (G.n_ind * 2);

	//Write records
	for (int l = 0 ; l < V.size() ; l ++) {

		//Get AC/AN
		int32_t vAC = 0, vAN = H.n_hap;
		for (int32_t h = 0 ; h < 2*G.n_ind ; h++) vAC += H.H_opt_var.get(l, h);

		//Variant information
		XW.writeInfo(V.vec_pos[l]->chr, V.vec_pos[l]->bp, V.vec_pos[l]->ref, V.vec_pos[l]->alt, V.vec_pos[l]->id, vAC, vAN);

		//Write haplotypes in BCF format
		if (hts_genotypes) {
			for (int32_t h = 0 ; h < 2*G.n_ind ; h++) output_buffer[h] = bcf_gt_phased(H.H_opt_var.get(l, h));
			XW.writeRecord(RECORD_BCFVCF_GENOTYPE, reinterpret_cast<char*>(output_buffer), 2 * G.n_ind * sizeof(int32_t));
		} else {
			for (int32_t h = 0 ; h < 2*G.n_ind ; h++) output_bitvector.set(h, H.H_opt_var.get(l, h));
			XW.writeRecord(RECORD_BINARY_HAPLOTYPE, reinterpret_cast<char*>(output_bitvector.bytes), output_bitvector.n_bytes);
		}

		//Verbose progress
		vrb.progress("  * VCF writing", (l+1)*1.0/V.size());
	}

	//Close
	free(output_buffer);
	XW.close();

	//Verbose writing
	if (hts_genotypes)
		vrb.bullet("VCF/BCF writing [N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)");
	else
		vrb.bullet("XCF writing [N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)");

	//Indexing
	vrb.bullet("Indexing files");
	if (bcf_index_build3(fname.c_str(), NULL, 14, nthreads) < 0) vrb.error("Fail to index file");
}
