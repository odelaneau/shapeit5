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
	if (!(bcf_sr_add_reader (sr, funphased.c_str()))) vrb.error("Problem opening index file for [" + funphased + "]");
	if (!(bcf_sr_add_reader (sr, fphased.c_str()))) vrb.error("Problem opening index file for [" + fphased + "]");

	//Sample processing // Needs to be improved to handle cases where sample lists do not properly overlap (in number and ordering)
	n_samples = bcf_hdr_nsamples(sr->readers[0].header);

	for (int i = 0 ; i < n_samples ; i ++)
		G.names.push_back(string(sr->readers[0].header->samples[i]));

	bcf1_t * line_phased, * line_unphased;
	int nset, vt = 0, vr = 0, vc = 0, vs = 0;
	int ngt_phased, *gt_arr_phased = NULL, ngt_arr_phased = 0;
	int ngt_unphased, *gt_arr_unphased = NULL, ngt_arr_unphased = 0;

	while (nset = bcf_sr_next_line (sr)) {
		line_unphased =  bcf_sr_get_line(sr, 0);
		line_phased =  bcf_sr_get_line(sr, 1);

		if (line_phased && line_phased->n_allele != 2) continue;
		if (line_unphased && line_phased->n_allele != 2) continue;

		char vartype = V.vec_full[vt]->type;
		bool minor =  V.vec_full[vt]->minor;

		if (line_phased) {
			assert(vartype == VARTYPE_SCAF);
			ngt_phased = bcf_get_genotypes(sr->readers[1].header, line_phased, &gt_arr_phased, &ngt_arr_phased); assert(ngt_phased == 2 * n_samples);
			for(int i = 0 ; i < 2 * n_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_phased[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_phased[i+1])==1);
				assert (gt_arr_phased[i+0] != bcf_gt_missing && gt_arr_phased[i+1] != bcf_gt_missing);
				H.Hvar.set(vs, i+0, a0);
				H.Hvar.set(vs, i+1, a1);
				n_scaffold_genotypes[a0 + a1] ++;
			}
		} else if (vartype == VARTYPE_COMM) {
			ngt_unphased = bcf_get_genotypes(sr->readers[0].header, line_unphased, &gt_arr_unphased, &ngt_arr_unphased); assert(ngt_unphased == 2 * n_samples);
			for(int i = 0 ; i < 2 * n_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_unphased[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_unphased[i+1])==1);
				bool mi = (gt_arr_unphased[i+0] == bcf_gt_missing || gt_arr_unphased[i+1] == bcf_gt_missing);
				if (mi) {
					G.setCommonMissing(vc, i);
					n_common_genotypes[3] ++;
				} else {
					G.setCommonGenotype(vc, i, a0+a1);
					n_common_genotypes[a0+a1] ++;
				}
			}
		} else {
			ngt_unphased = bcf_get_genotypes(sr->readers[0].header, line_unphased, &gt_arr_unphased, &ngt_arr_unphased); assert(ngt_unphased == 2 * n_samples);
			for(int i = 0 ; i < 2 * n_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_unphased[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_unphased[i+1])==1);
				bool mi = (gt_arr_unphased[i+0] == bcf_gt_missing || gt_arr_unphased[i+1] == bcf_gt_missing);
				if (mi) {
					G.pushRareMissing(vr, i);
					n_rare_genotypes[3] ++;
				} else if ((a0+a1) == 1) {
					G.pushRareHet(vr, i);
					n_rare_genotypes[1] ++;
				} else if (a0 == minor) {
					G.pushRareHom(vr, i);
					n_rare_genotypes[a0*2] ++;
				} else n_rare_genotypes[a0*2] ++;
			}
		}

		switch (vartype) {
		case VARTYPE_SCAF :	vs ++; vt ++; break;
		case VARTYPE_COMM :	vc ++; vt ++; break;
		case VARTYPE_RARE :	vr ++; vt ++; break;
		default: vrb.error("Unrecognized variant type!");
		}
		vrb.progress("  * VCF/BCF parsing", vt * 1.0 / n_total_variants);
	}
	free(gt_arr_unphased);
	free(gt_arr_phased);
	bcf_sr_destroy(sr);

	// Report
	vrb.bullet("VCF/BCF parsing ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("  + Scaffold : R/R=" + stb.str(n_scaffold_genotypes[0]) + " R/A=" + stb.str(n_scaffold_genotypes[1]) + " A/A=" + stb.str(n_scaffold_genotypes[2]) + " ./.=" + stb.str(n_scaffold_genotypes[3]));
	vrb.bullet("  + Common   : R/R=" + stb.str(n_common_genotypes[0]) + " R/A=" + stb.str(n_common_genotypes[1]) + " A/A=" + stb.str(n_common_genotypes[2]) + " ./.=" + stb.str(n_common_genotypes[3]));
	vrb.bullet("  + Rare     : R/R=" + stb.str(n_rare_genotypes[0]) + " R/A=" + stb.str(n_rare_genotypes[1]) + " A/A=" + stb.str(n_rare_genotypes[2]) + " ./.=" + stb.str(n_rare_genotypes[3]));
}
