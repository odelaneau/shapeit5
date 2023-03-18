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

void genotype_reader::readGenotypes() {
	tac.clock();
	vrb.wait("  * VCF/BCF parsing");

	//Initialize VCF/BCF reader(s)
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "]");

	//Opening file(s)
	for (int f = 0 ; f < 3 ; f ++) if (panels[f] && !(bcf_sr_add_reader (sr, filenames[f].c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + filenames[f] + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + filenames[f] + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + filenames[f] + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + filenames[f] + "]: unknown error"); break;
		}
	}

	//Main sample IDs processing
	std::map < std::string, int > map_names;
	for (int i = 0 ; i < n_main_samples ; i ++) {
		G.vecG[i]->name = std::string(sr->readers[0].header->samples[i]);
		map_names.insert(std::pair < std::string, int > (G.vecG[i]->name, i));
	}

	//Scaffold sample IDs processing
	int n_scaf_samples = 0;
	std::vector < int > mappingS2M;
	if (panels[2]) {
		int n_with_scaffold = 0;
		n_scaf_samples = bcf_hdr_nsamples(sr->readers[panels[1]+panels[2]].header);
		mappingS2M = std::vector < int > (n_scaf_samples, -1);
		std::vector < std::string > tokens_tmp;
		for (int i = 0 ; i < n_scaf_samples ; i ++) {
			std::string buffer_tmp = std::string(sr->readers[panels[1]+panels[2]].header->samples[i]);
			stb.split(buffer_tmp, tokens_tmp, "_");
			std::string scaf_name = tokens_tmp[0];
			std::map < std::string, int > :: iterator it = map_names.find(scaf_name);
			if (it != map_names.end()) {
				mappingS2M[i] = it->second;
				n_with_scaffold++;
			}
		}
		vrb.bullet(stb.str(n_with_scaffold) + " samples with scaffold");
	}

	//Parsing VCF/BCF
	unsigned int i_variant_total = 0, i_variant_kept = 0, nset = 0;
	int ngt_main, *gt_arr_main = NULL, ngt_arr_main = 0;
	int ngt_scaf, *gt_arr_scaf = NULL, ngt_arr_scaf = 0;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_main, * line_ref, * line_scaf;
	while ((nset = bcf_sr_next_line (sr))) {
		if (variant_mask[i_variant_total]) {

			//Retrieve data in main file
			line_main = bcf_sr_get_line(sr, 0);
			ngt_main = bcf_get_genotypes(sr->readers[0].header, line_main, &gt_arr_main, &ngt_arr_main); assert(ngt_main == 2 * n_main_samples);

			//Retrieve data in reference file if necessary
			if (panels[1]) {
				line_ref = bcf_sr_get_line(sr, 1);
				ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref); assert(ngt_ref == 2 * n_ref_samples);
			}

			//Retrieve data in scaffold file if necessary and available
			if (panels[2]) {
				line_scaf = bcf_sr_get_line(sr, panels[1]+panels[2]);
				if (line_scaf) {
					ngt_scaf = bcf_get_genotypes(sr->readers[panels[1]+panels[2]].header, line_scaf, &gt_arr_scaf, &ngt_arr_scaf);
					assert(ngt_scaf == 2 * n_scaf_samples);
				}
			}

			//Process main data
			for(int i = 0 ; i < 2 * n_main_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_main[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_main[i+1])==1);
				bool mi = (gt_arr_main[i+0] == bcf_gt_missing || gt_arr_main[i+1] == bcf_gt_missing);
				bool he = !mi && a0 != a1;
				bool ho = !mi && a0 == a1;
				if (a0) VAR_SET_HAP0(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
				if (a1) VAR_SET_HAP1(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
				if (mi) VAR_SET_MIS(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
				if (he) VAR_SET_HET(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
				if (mi) {
					V.vec_pos[i_variant_kept]->cmis++;
					n_genotypes[3] ++;
				} else {
					V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1);
					V.vec_pos[i_variant_kept]->calt += a0+a1;
					n_genotypes[a0+a1] ++;
				}
			}

			//Process reference data
			if (panels[1]) {
				for(int i = 0 ; i < 2 * n_ref_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(gt_arr_ref[i+0])==1);
					bool a1 = (bcf_gt_allele(gt_arr_ref[i+1])==1);
					if (gt_arr_ref[i+0] == bcf_gt_missing || gt_arr_ref[i+1] == bcf_gt_missing) vrb.error("Missing genotype(s) in reference panel");
					if (!bcf_gt_is_phased(gt_arr_ref[i+1])) vrb.error("Unphased genotype(s) in reference panel");
					H.H_opt_hap.set(i+2*n_main_samples+0, i_variant_kept, a0);
					H.H_opt_hap.set(i+2*n_main_samples+1, i_variant_kept, a1);
					V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1);
					V.vec_pos[i_variant_kept]->calt += a0+a1;
				}
			}

			//Process scaffold data
			if (panels[2] && line_scaf) {
				for(int i = 0 ; i < 2 * n_scaf_samples ; i += 2) {
					int ind = mappingS2M[DIV2(i)];
					if (ind >= 0) {
						bool sa0 = (bcf_gt_allele(gt_arr_scaf[i+0])==1);
						bool sa1 = (bcf_gt_allele(gt_arr_scaf[i+1])==1);
						bool sph = (bcf_gt_is_phased(gt_arr_scaf[i+0]) || bcf_gt_is_phased(gt_arr_scaf[i+1]));
						bool smi = (gt_arr_scaf[i+0] == bcf_gt_missing || gt_arr_scaf[i+1] == bcf_gt_missing);
						if ((sa0 != sa1) && !smi && sph && VAR_GET_HET(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)])) {
							VAR_SET_SCA(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
							sa0?VAR_SET_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
							sa1?VAR_SET_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
							n_genotypes[4] ++;
						}
					}
				}
			}
			vrb.progress("  * VCF/BCF parsing", i_variant_kept * 1.0 / n_variants);
			i_variant_kept ++;
		}
		i_variant_total++;
	}
	free(gt_arr_main);
	if (panels[1]) free(gt_arr_ref);
	if (panels[2]) free(gt_arr_scaf);
	bcf_sr_destroy(sr);

	// Report
	unsigned long n_genotypes_total = accumulate(n_genotypes.begin(), n_genotypes.begin() + 4, 0UL);
	std::string str0 = "0/0=" + stb.str(n_genotypes[0]*100.0/n_genotypes_total, 3) + "%";
	std::string str1 = "0/1=" + stb.str(n_genotypes[1]*100.0/n_genotypes_total, 3) + "%";
	std::string str2 = "1/1=" + stb.str(n_genotypes[2]*100.0/n_genotypes_total, 3) + "%";
	std::string str3 = "./.=" + stb.str(n_genotypes[3]*100.0/n_genotypes_total, 3) + "%";
	std::string str4 = "0|1=" + stb.str(n_genotypes[4]*100.0/n_genotypes_total, 3) + "%";
	std::string str5 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing ["+str0+" / "+str1+" / "+str2+" / "+str3+" / "+str4+"] ("+str5+")");
}
