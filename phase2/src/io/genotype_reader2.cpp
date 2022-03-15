////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <io/genotype_reader.h>

//**********************************************************************************//
//								ONE VCF/BCF PROCESSED								//
//								1. main genotype data								//
//**********************************************************************************//
void genotype_reader::readGenotypes0(string funphased) {
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader(sr, funphased.c_str());
	for (int i = 0 ; i < n_target_samples ; i ++) G.vecG[i]->name = string(sr->readers[0].header->samples[i]);
	bcf1_t * line;
	int ngt_target, *gt_arr_target = NULL, ngt_arr_target = 0;
	int nps_target, *ps_arr_target = NULL, nps_arr_target = 0;

	unsigned int i_common_variant = 0, i_total_variant = 0;

	while(bcf_sr_next_line (sr)) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {

			bool minor_allele = V.vec_full[i_total_variant]->minor;
			ngt_target = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr_target, &ngt_arr_target); assert(ngt_target == 2 * n_target_samples);
			if (use_PS_field) {
				nps_target = bcf_get_format_int32(sr->readers[0].header, line, "PS", &ps_arr_target, &nps_arr_target);
				setPScodes(ps_arr_target, nps_target);
			}

			for(int i = 0 ; i < 2 * n_target_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_target[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_target[i+1])==1);
				bool mi = (gt_arr_target[i+0] == bcf_gt_missing || gt_arr_target[i+1] == bcf_gt_missing);
				bool he = !mi && a0 != a1;
				bool ho = !mi && a0 == a1;
				bool ps = he && use_PS_field;
				bool ph = (bcf_gt_is_phased(gt_arr_target[i+0]) || bcf_gt_is_phased(gt_arr_target[i+1])) && he && PScodes.size() > 0;

				if (flag_common[i_total_variant]) {
					if (a0) VAR_SET_HAP0(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (a1) VAR_SET_HAP1(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (mi) VAR_SET_MIS(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (he) VAR_SET_HET(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (ps) G.vecG[DIV2(i)]->pushPS(a0, a1, ph?PScodes[i/2]:0);
				} else if (mi || (!mi && ((a0 == minor_allele) || (a1 == minor_allele)))) {
					unsigned int rare_idx = G.vecG[DIV2(i)]->RareIndexes.size();
					if (!MOD2(rare_idx)) G.vecG[DIV2(i)]->RareVariants.push_back(0);
					if (a0) VAR_SET_HAP0(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					if (a1) VAR_SET_HAP1(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					if (mi) VAR_SET_MIS(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					if (he) VAR_SET_HET(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					G.vecG[DIV2(i)]->RareIndexes.push_back(i_total_variant);
				}

				if (!mi) {
					a0?V.vec_full[i_total_variant]->calt++:V.vec_full[i_total_variant]->cref++;
					a1?V.vec_full[i_total_variant]->calt++:V.vec_full[i_total_variant]->cref++;
				} else V.vec_full[i_total_variant]->cmis++;

				n_genotypes[0] ++;
				n_genotypes[1] += ho && !a0;
				n_genotypes[2] += he;
				n_genotypes[3] += ho && a0;
				n_genotypes[4] += mi;
				n_genotypes[5] += ph && ps;
			}
			i_common_variant += flag_common[i_total_variant];
			i_total_variant ++;
			vrb.progress("  * VCF/BCF parsing", i_total_variant*1.0/n_total_variants);
		}
	}
	free(gt_arr_target);
	if (ps_arr_target) free(ps_arr_target);
	bcf_sr_destroy(sr);
	// Report
	string str0 = "RR=" + stb.str(n_genotypes[1]*100.0/n_genotypes[0], 1) + "%";
	string str1 = "RA=" + stb.str(n_genotypes[2]*100.0/n_genotypes[0], 1) + "%" + (use_PS_field?(" / PS=" + stb.str(n_genotypes[5]*100.0/n_genotypes[0], 3) + "%"):(""));
	string str2 = "AA=" + stb.str(n_genotypes[3]*100.0/n_genotypes[0], 1) + "%";
	string str3 = "MI=" + stb.str(n_genotypes[4]*100.0/n_genotypes[0], 1) + "%";
	string strt = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing ["+str0+" / "+str1+" / "+str2+" / "+str3+"] ("+strt+")");
}

//**********************************************************************************//
//								TWO VCF/BCF PROCESSED								//
//								1. main genotype data								//
//								2. scaffold haplotype data							//
//**********************************************************************************//
void genotype_reader::readGenotypes1(string funphased, string fphased) {
	tac.clock();

	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader(sr, funphased.c_str());
	if (!bcf_sr_add_reader (sr, fphased.c_str())) vrb.error("Problem opening index file for [" + fphased + "]");

	// Mapping scaffolded samples
	map < string, int > map_names;
	for (int i = 0 ; i < n_target_samples ; i ++) {
		G.vecG[i]->name = string(sr->readers[0].header->samples[i]);
		map_names.insert(pair < string, int > (G.vecG[i]->name, i));
	}
	int n_scaffold_samples = bcf_hdr_nsamples(sr->readers[1].header);
	vector < int > mappingS2T = vector < int > (n_scaffold_samples, -1);
	for (int i = 0 ; i < n_scaffold_samples ; i ++) {
		string scaffold_name = string(sr->readers[1].header->samples[i]);
		map < string, int > :: iterator it = map_names.find(scaffold_name);
		if (it != map_names.end()) mappingS2T[i] = it->second;
	}

	bcf1_t * line_target, * line_scaffold;
	int ngt_target, *gt_arr_target = NULL, ngt_arr_target = 0;
	int nps_target, *ps_arr_target = NULL, nps_arr_target = 0;
	int ngt_scaffo, *gt_arr_scaffo = NULL, ngt_arr_scaffo = 0;

	unsigned int i_common_variant = 0, i_total_variant = 0;

	while (bcf_sr_next_line (sr)) {
		line_target =  bcf_sr_get_line(sr, 0);
		if (line_target->n_allele == 2) {

			//Fetch data in target VCF/BCF file
			bool minor_allele = V.vec_full[i_total_variant]->minor;
			ngt_target = bcf_get_genotypes(sr->readers[0].header, line_target, &gt_arr_target, &ngt_arr_target); assert(ngt_target == 2 * n_target_samples);
			if (use_PS_field) {
				nps_target = bcf_get_format_int32(sr->readers[0].header, line_target, "PS", &ps_arr_target, &nps_arr_target);
				setPScodes(ps_arr_target, nps_target);
			}

			for(int i = 0 ; i < 2 * n_target_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_target[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_target[i+1])==1);
				bool mi = (gt_arr_target[i+0] == bcf_gt_missing || gt_arr_target[i+1] == bcf_gt_missing);
				bool he = !mi && a0 != a1;
				bool ho = !mi && a0 == a1;
				bool ps = he && use_PS_field;
				bool ph = (bcf_gt_is_phased(gt_arr_target[i+0]) || bcf_gt_is_phased(gt_arr_target[i+1])) && he && PScodes.size() > 0;

				if (flag_common[i_total_variant]) {
					if (a0) VAR_SET_HAP0(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (a1) VAR_SET_HAP1(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (mi) VAR_SET_MIS(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (he) VAR_SET_HET(MOD2(i_common_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_common_variant)]);
					if (ps) G.vecG[DIV2(i)]->pushPS(a0, a1, ph?PScodes[i/2]:0);
				} else if (mi || (!mi && ((a0 == minor_allele) || (a1 == minor_allele)))) {
					unsigned int rare_idx = G.vecG[DIV2(i)]->RareIndexes.size();
					if (!MOD2(rare_idx)) G.vecG[DIV2(i)]->RareVariants.push_back(0);
					if (a0) VAR_SET_HAP0(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					if (a1) VAR_SET_HAP1(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					if (mi) VAR_SET_MIS(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					if (he) VAR_SET_HET(MOD2(rare_idx), G.vecG[DIV2(i)]->RareVariants[DIV2(rare_idx)]);
					G.vecG[DIV2(i)]->RareIndexes.push_back(i_total_variant);
				}

				if (!mi) {
					a0?V.vec_full[i_total_variant]->calt++:V.vec_full[i_total_variant]->cref++;
					a1?V.vec_full[i_total_variant]->calt++:V.vec_full[i_total_variant]->cref++;
				} else V.vec_full[i_total_variant]->cmis++;

				n_genotypes[0] ++;
				n_genotypes[1] += ho && !a0;
				n_genotypes[2] += he;
				n_genotypes[3] += ho && a0;
				n_genotypes[4] += mi;
				n_genotypes[5] += ph && ps;
			}

			//Fetch data in scaffold VCF/BCF file
			if (flag_common[i_total_variant] && (line_scaffold = bcf_sr_get_line(sr, 1))) {
				ngt_scaffo = bcf_get_genotypes(sr->readers[1].header, line_scaffold, &gt_arr_scaffo, &ngt_arr_scaffo); assert(ngt_scaffo == 2 * n_scaffold_samples);
				for(int i = 0 ; i < 2 * n_scaffold_samples ; i += 2) {
					int index_ind = mappingS2T[DIV2(i)];
					if (index_ind >= 0) {
						bool s0 = (bcf_gt_allele(gt_arr_scaffo[i+0])==1);
						bool s1 = (bcf_gt_allele(gt_arr_scaffo[i+1])==1);
						bool he = (s0 != s1);
						bool ph = (bcf_gt_is_phased(gt_arr_scaffo[i+0]) || bcf_gt_is_phased(gt_arr_scaffo[i+1]));
						bool mi = (gt_arr_scaffo[i+0] == bcf_gt_missing || gt_arr_scaffo[i+1] == bcf_gt_missing);
						if (he && !mi && ph) {
							bool a0 = VAR_GET_HAP0(MOD2(i_common_variant), G.vecG[index_ind]->Variants[DIV2(i_common_variant)]);
							bool a1 = VAR_GET_HAP1(MOD2(i_common_variant), G.vecG[index_ind]->Variants[DIV2(i_common_variant)]);
							if (a0!=a1) {	//Only scaffold common variants / target and scaffold genotypes need to be het and match
								VAR_SET_SCA(MOD2(i_common_variant), G.vecG[index_ind]->Variants[DIV2(i_common_variant)]);
								s0?VAR_SET_HAP0(MOD2(i_common_variant), G.vecG[index_ind]->Variants[DIV2(i_common_variant)]):VAR_CLR_HAP0(MOD2(i_common_variant), G.vecG[index_ind]->Variants[DIV2(i_common_variant)]);
								s1?VAR_SET_HAP1(MOD2(i_common_variant), G.vecG[index_ind]->Variants[DIV2(i_common_variant)]):VAR_CLR_HAP1(MOD2(i_common_variant), G.vecG[index_ind]->Variants[DIV2(i_common_variant)]);
								n_genotypes[6] ++;
							}
						}
					}
				}
			}

			i_common_variant += flag_common[i_total_variant];
			i_total_variant ++;
			vrb.progress("  * VCF/BCF parsing", i_total_variant*1.0/n_total_variants);
		}
	}
	free(gt_arr_target);
	if (ps_arr_target) free(ps_arr_target);
	bcf_sr_destroy(sr);
	// Report
	string str0 = "RR=" + stb.str(n_genotypes[1]*100.0/n_genotypes[0], 1) + "%";
	string str1 = "RA=" + stb.str(n_genotypes[2]*100.0/n_genotypes[0], 1) + "%" + (use_PS_field?(" / PS=" + stb.str(n_genotypes[5]*100.0/n_genotypes[0], 3) + "%"):(""));
	string str2 = "AA=" + stb.str(n_genotypes[3]*100.0/n_genotypes[0], 1) + "%";
	string str3 = "MI=" + stb.str(n_genotypes[4]*100.0/n_genotypes[0], 1) + "%";
	string str4 = "SC=" + stb.str(n_genotypes[6]*100.0/n_genotypes[0], 1) + "%";
	string strt = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing ["+str0+" / "+str1+" / "+str2+" / "+str3+" / "+str4+"] ("+strt+")");
}
