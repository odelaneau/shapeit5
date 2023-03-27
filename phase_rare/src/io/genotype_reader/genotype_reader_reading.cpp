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
#include <containers/bitvector.h>
#include <io/genotype_reader/genotype_reader_header.h>

void genotype_reader::readGenotypes() {
	tac.clock();
	vrb.wait("  * BCF parsing");

	//Initialize XCF reader
	xcf_reader XR (scaffold_region, nthreads);

	//Open files
	uint32_t idx_unphased = XR.addFile(funphased);
	uint32_t idx_scaffold = XR.addFile(fphased);

	//Allocate buffers
	int32_t * input_buffer = (int32_t*)malloc(2 * n_samples * sizeof(int32_t));
	bitvector input_bitvector (2 * n_samples);

	//Parse file
	int32_t vt = 0, vr = 0, vs = 0;
	while (XR.nextRecord()) {

		//CASE1: Scaffold variant ...
		if (XR.hasRecord(idx_scaffold)) {

			//Checks
			assert(V.vec_full[vt]->type == VARTYPE_SCAF);
			int32_t scaf_type = XR.typeRecord(idx_scaffold);

			//... in BCF format
			if (scaf_type == RECORD_BCFVCF_GENOTYPE) {
				XR.readRecord(idx_scaffold, reinterpret_cast< char** > (&input_buffer));
				for(int i = 0 ; i < 2 * n_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(input_buffer[i+0])==1);
					bool a1 = (bcf_gt_allele(input_buffer[i+1])==1);
					if (input_buffer[i+0] == bcf_gt_missing || input_buffer[i+1] == bcf_gt_missing) vrb.error("Missing data is not allowed in scaffold");
					V.vec_full[vt]->cref += 2-(a0+a1);
					V.vec_full[vt]->calt += a0+a1;
					H.Hvar.set(vs, i+0, a0);
					H.Hvar.set(vs, i+1, a1);
					n_scaffold_genotypes[a0 + a1] ++;
				}
				vs++; vt ++;
			}

			//... in Binary format
			else if (scaf_type == RECORD_BINARY_HAPLOTYPE) {
				XR.readRecord(idx_scaffold, reinterpret_cast< char** > (&input_bitvector.bytes));
				for(int i = 0 ; i < 2 * n_samples ; i += 2) {
					bool a0 = input_bitvector.get(i+0);
					bool a1 = input_bitvector.get(i+1);
					V.vec_full[vt]->cref += 2-(a0+a1);
					V.vec_full[vt]->calt += a0+a1;
					H.Hvar.set(vs, i+0, a0);
					H.Hvar.set(vs, i+1, a1);
					n_scaffold_genotypes[a0 + a1] ++;
				}
				vs++; vt ++;
			}

			//... unsupported format for scaffold
			else vrb.error("Record format [" + stb.str(scaf_type) + "] is not supported for scaffold data");
		}

		//CASE2: Rare variant ...
		else if (XR.hasRecord(idx_unphased) && XR.pos >= input_start && XR.pos <= input_stop) {

			//Checks
			assert(V.vec_full[vt]->type == VARTYPE_RARE);
			int32_t unph_type = XR.typeRecord(idx_unphased);
			bool minor = (XR.getAF(idx_unphased) < 0.5f);

			// ... in BCF format
			if (unph_type == RECORD_BCFVCF_GENOTYPE) {
				XR.readRecord(idx_unphased, reinterpret_cast< char** > (&input_buffer));
				for(int i = 0 ; i < 2 * n_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(input_buffer[i+0])==1);
					bool a1 = (bcf_gt_allele(input_buffer[i+1])==1);
					bool mi = (input_buffer[i+0] == bcf_gt_missing || input_buffer[i+1] == bcf_gt_missing);
					if (mi) V.vec_full[vt]->cmis ++;
					else { V.vec_full[vt]->cref += 2-(a0+a1); V.vec_full[vt]->calt += a0+a1; }

					if (mi) {
						G.pushRareMissing(vr, i/2, !minor);
						n_rare_genotypes[3] ++;
					} else if ((a0+a1) == 1) {
						G.pushRareUnphasedHet(vr, i/2);
						n_rare_genotypes[1] ++;
					} else if (a0 == minor) {
						G.pushRareHom(vr, i/2, !minor);
						n_rare_genotypes[a0*2] ++;
					} else n_rare_genotypes[a1*2] ++;
				}
				vr++; vt ++;
			}

			// ... in Binary format
			else if (unph_type == RECORD_BINARY_GENOTYPE) {
				XR.readRecord(idx_unphased, reinterpret_cast< char** > (&input_bitvector.bytes));
				for(int i = 0 ; i < 2 * n_samples ; i += 2) {
					bool a0 = input_bitvector.get(i+0);
					bool a1 = input_bitvector.get(i+1);
					bool mi = a0 && !a1;
					if (mi) V.vec_full[vt]->cmis ++;
					else { V.vec_full[vt]->cref += 2-(a0+a1); V.vec_full[vt]->calt += a0+a1; }
					if (mi) {
						G.pushRareMissing(vr, i/2, !minor);
						n_rare_genotypes[3] ++;
					} else if ((a0+a1) == 1) {
						G.pushRareUnphasedHet(vr, i/2);
						n_rare_genotypes[1] ++;
					} else if (a0 == minor) {
						G.pushRareHom(vr, i/2, !minor);
						n_rare_genotypes[a0*2] ++;
					} else n_rare_genotypes[a1*2] ++;
				}
				vr++; vt ++;
			}

			// ... in Sparse format
			else if (unph_type == RECORD_SPARSE_GENOTYPE) {
				//Standard version, involves decoding/encoding of sparse data. Could be improve by direct copy, but not sure it would drastically speed up things 				//G.GRvar_genotypes[vr].reserve(XR.sizeRecord(idx_unphased) / sizeof(int32_t));		//XR.readRecord(idx_unphased, reinterpret_cast< char** > (&G.GRvar_genotypes[vr].data()));
				int32_t n_elements = XR.readRecord(idx_unphased, reinterpret_cast< char** > (&input_buffer)) / sizeof(int32_t);
				G.GRvar_genotypes[vr].reserve(n_elements);
				for(int r = 0 ; r < n_elements ; r++) n_rare_genotypes[G.pushRare(vr, input_buffer[r])] ++;
				vr++; vt ++;
			}

			//... unsupported format for scaffold
			else vrb.error("Record format [" + stb.str(unph_type) + "] is not supported for unphased data");
		}

		//Verbose
		vrb.progress("  * BCF parsing", vt * 1.0f / n_total_variants);
	}

	//Close
	free(input_buffer);
	XR.close();

	// Report
	vrb.bullet("Plain VCF/BCF parsing ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("Scaffold : [0/0=" + stb.str(n_scaffold_genotypes[0]) + " 0/1=" + stb.str(n_scaffold_genotypes[1]) + " 1/1=" + stb.str(n_scaffold_genotypes[2]) + "]");
	vrb.bullet2("Rare     : [0/0=" + stb.str(n_rare_genotypes[0]) + " 0/1=" + stb.str(n_rare_genotypes[1]) + " 1/1=" + stb.str(n_rare_genotypes[2]) + " ./.=" + stb.str(n_rare_genotypes[3]) + "]");
}
