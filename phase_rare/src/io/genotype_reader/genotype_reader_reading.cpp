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

void genotype_reader::readGenotypesPlain() {
	tac.clock();
	vrb.wait("  * Plain VCF/BCF parsing");

	//Initialize VCF/BCF reader(s)
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, scaffold_region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + scaffold_region + "]");

	//Opening file(s)
	if (!(bcf_sr_add_reader (sr, funphased.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + funphased + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + funphased + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + funphased + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + funphased + "]: unknown error"); break;
		}
	}
	if (!(bcf_sr_add_reader (sr, fphased.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + fphased + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + fphased + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + fphased + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + fphased + "]: unknown error"); break;
		}
	}

	//Sample processing
	n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples ; i ++) G.names.push_back(std::string(sr->readers[0].header->samples[i]));

	//Check if sample lists match
	if (bcf_hdr_nsamples(sr->readers[0].header) != bcf_hdr_nsamples(sr->readers[1].header)) vrb.error("Different number of samples in the specified BCF/BCF files");
	for (int i = 0 ; i < n_samples ; i ++) if (std::string(sr->readers[1].header->samples[i]) != G.names[i]) vrb.error("Different sample IDs in the specified BCF/BCF files");

	//Parse file
	bcf1_t * line_phased = NULL, * line_unphased = NULL;
	int nset, vt = 0, vr = 0, vs = 0;
	int ngt_phased, *gt_arr_phased = NULL, ngt_arr_phased = 0;
	int ngt_unphased, *gt_arr_unphased = NULL, ngt_arr_unphased = 0;
	while ((nset = bcf_sr_next_line (sr))) {
		line_unphased =  bcf_sr_get_line(sr, 0);
		line_phased =  bcf_sr_get_line(sr, 1);

		if (line_phased && line_phased->n_allele != 2) continue;
		if (line_unphased && line_unphased->n_allele != 2) continue;

		if (line_phased) {
			assert(V.vec_full[vt]->type == VARTYPE_SCAF);
			ngt_phased = bcf_get_genotypes(sr->readers[1].header, line_phased, &gt_arr_phased, &ngt_arr_phased); assert(ngt_phased == 2 * n_samples);
			for(int i = 0 ; i < 2 * n_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_phased[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_phased[i+1])==1);
				assert (gt_arr_phased[i+0] != bcf_gt_missing && gt_arr_phased[i+1] != bcf_gt_missing);
				V.vec_full[vt]->cref += 2-(a0+a1);
				V.vec_full[vt]->calt += a0+a1;
				H.Hvar.set(vs, i+0, a0);
				H.Hvar.set(vs, i+1, a1);
				n_scaffold_genotypes[a0 + a1] ++;
			}
			vs++; vt ++;
		} else {
			int pos = line_unphased->pos + 1;
			if (pos >= input_start && pos <= input_stop) {
				if (V.vec_full[vt]->type == VARTYPE_RARE) {
					bool minor =  V.vec_full[vt]->minor;
					ngt_unphased = bcf_get_genotypes(sr->readers[0].header, line_unphased, &gt_arr_unphased, &ngt_arr_unphased); assert(ngt_unphased == 2 * n_samples);
					for(int i = 0 ; i < 2 * n_samples ; i += 2) {
						bool a0 = (bcf_gt_allele(gt_arr_unphased[i+0])==1);
						bool a1 = (bcf_gt_allele(gt_arr_unphased[i+1])==1);
						bool mi = (gt_arr_unphased[i+0] == bcf_gt_missing || gt_arr_unphased[i+1] == bcf_gt_missing);
						if (mi) V.vec_full[vt]->cmis ++;
						else {
							V.vec_full[vt]->cref += 2-(a0+a1);
							V.vec_full[vt]->calt += a0+a1;
						}

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
			}
		}
		vrb.progress("  * Plain VCF/BCF parsing", vt * 1.0f / n_total_variants);
	}


	free(gt_arr_unphased);
	free(gt_arr_phased);
	bcf_sr_destroy(sr);

	// Report
	vrb.bullet("Plain VCF/BCF parsing ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("Scaffold : R/R=" + stb.str(n_scaffold_genotypes[0]) + " R/A=" + stb.str(n_scaffold_genotypes[1]) + " A/A=" + stb.str(n_scaffold_genotypes[2]) + " ./.=" + stb.str(n_scaffold_genotypes[3]));
	vrb.bullet2("Rare     : R/R=" + stb.str(n_rare_genotypes[0]) + " R/A=" + stb.str(n_rare_genotypes[1]) + " A/A=" + stb.str(n_rare_genotypes[2]) + " ./.=" + stb.str(n_rare_genotypes[3]));
}


void genotype_reader::readGenotypesSparse() {
	tac.clock();
	vrb.wait("  * Sparse VCF/BCF parsing");

	//Initialize VCF/BCF reader(s)
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, scaffold_region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + scaffold_region + "]");

	//Opening file(s)
	if (!(bcf_sr_add_reader (sr, fphased.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + fphased + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + fphased + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + fphased + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + fphased + "]: unknown error"); break;
		}
	}
	if (!(bcf_sr_add_reader (sr, funphased.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + funphased + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + funphased + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + funphased + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + funphased + "]: unknown error"); break;
		}
	}

	//Opening sparse BIN file
	int id = bcf_hdr_id2int(sr->readers[1].header, BCF_DT_ID, "SGEN");
	if (!bcf_hdr_idinfo_exists(sr->readers[1].header, BCF_HL_INFO, id)) vrb.error("Cannot retrieve SGEN flag in VCF/BCF for sparse rare variants.");
	std::string fbinary = stb.get_name_from_vcf(funphased) + ".bin";
	std::ifstream fp_binary (fbinary, std::ios::in | std::ios::binary);
	if (!fp_binary) vrb.error("Cannot open " + fbinary + " for reading, check permissions");

	//Sample processing
	n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples ; i ++) G.names.push_back(std::string(sr->readers[0].header->samples[i]));

	bcf1_t * line_phased = NULL, * line_unphased = NULL;
	int nset, vt = 0, vr = 0, vs = 0;
	int nsk = 0, rsk = 0, *vsk = NULL;
	int ngt_phased, *gt_arr_phased = NULL, ngt_arr_phased = 0;
	unsigned int * rg_buffer = (unsigned int *)malloc(n_samples * sizeof(unsigned int));

	while ((nset = bcf_sr_next_line (sr))) {
		line_phased =  bcf_sr_get_line(sr, 0);
		line_unphased =  bcf_sr_get_line(sr, 1);

		//Skip non bi-allelic
		if (line_phased && line_phased->n_allele != 2) continue;
		if (line_unphased && line_unphased->n_allele != 2) continue;

		//Check for duplicate variants in the two file
		if (nset == 2) vrb.error("Duplicate variant at pos= " + stb.str(line_phased->pos+1));

		if (line_phased) {
			assert(V.vec_full[vt]->type == VARTYPE_SCAF);
			ngt_phased = bcf_get_genotypes(sr->readers[0].header, line_phased, &gt_arr_phased, &ngt_arr_phased); assert(ngt_phased == 2 * n_samples);
			for(int i = 0 ; i < 2 * n_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_phased[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_phased[i+1])==1);
				assert (gt_arr_phased[i+0] != bcf_gt_missing && gt_arr_phased[i+1] != bcf_gt_missing);
				V.vec_full[vt]->cref += 2-(a0+a1);
				V.vec_full[vt]->calt += a0+a1;
				H.Hvar.set(vs, i+0, a0);
				H.Hvar.set(vs, i+1, a1);
				n_scaffold_genotypes[a0 + a1] ++;
			}
			vs++; vt ++;
		} else {
			int pos = line_unphased->pos + 1;
			if (pos >= input_start && pos <= input_stop) {
				assert(V.vec_full[vt]->type == VARTYPE_RARE);
				bool minor =  V.vec_full[vt]->minor;

				rsk = bcf_get_info_int32(sr->readers[1].header, line_unphased, "SGEN", &vsk, &nsk); if (nsk!=3) vrb.error("SGEN field is needed in rare file");

				uint64_t sgenotypes = vsk[0];
				sgenotypes *= MOD30BITS;
				sgenotypes += vsk[1];
				uint32_t ngenotypes = vsk[2];

				fp_binary.seekg(sgenotypes*sizeof(unsigned int), fp_binary.beg);
				fp_binary.read((char *)rg_buffer, ngenotypes * sizeof(unsigned int));

				for (int r = 0 ; r < ngenotypes ; r++) n_rare_genotypes[G.pushRare(vr, rg_buffer[r])] ++;
				n_rare_genotypes[2*(1-minor)] += n_samples-ngenotypes;

				vr++; vt ++;
			}
		}
		vrb.progress("  * Sparse VCF/BCF parsing", vt * 1.0 / n_total_variants);
	}
	free(rg_buffer);
	free(gt_arr_phased);
	free(vsk);
	bcf_sr_destroy(sr);
	fp_binary.close();

	// Report
	vrb.bullet("Sparse VCF/BCF parsing ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("Scaffold : R/R=" + stb.str(n_scaffold_genotypes[0]) + " R/A=" + stb.str(n_scaffold_genotypes[1]) + " A/A=" + stb.str(n_scaffold_genotypes[2]) + " ./.=" + stb.str(n_scaffold_genotypes[3]));
	vrb.bullet2("Rare     : R/R=" + stb.str(n_rare_genotypes[0]) + " R/A=" + stb.str(n_rare_genotypes[1]) + " A/A=" + stb.str(n_rare_genotypes[2]) + " ./.=" + stb.str(n_rare_genotypes[3]));
}
