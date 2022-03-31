/*******************************************************************************
 * Copyright (C) 2018-2022 Olivier Delaneau, University of Lausanne
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

void genotype_reader::scanGenotypes() {
	tac.clock();
	vrb.wait("  * VCF/BCF scanning");

	//Initialize VCF/BCF reader(s)
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "]");

	//Opening file(s)
	for (int f = 0 ; f < 2 ; f ++) if (panels[f] && !(bcf_sr_add_reader (sr, filenames[f].c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + filenames[f] + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + filenames[f] + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + filenames[f] + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + filenames[f] + "]: unknown error"); break;
		}
	}

	//Sample processing
	n_main_samples = bcf_hdr_nsamples(sr->readers[0].header);
	n_ref_samples = panels[1] ? bcf_hdr_nsamples(sr->readers[1].header) : 0;
	sample_mask = vector < bool > (n_main_samples + n_ref_samples, true);	//SET THIS UP!

	bcf1_t * line_main;
	int nset, n_variants_noverlap = 0, n_variants_multi = 0, n_variants_notsnp = 0, n_variants_rare = 0;
	int rAC_main=0, nAC_main=0, *vAC_main=NULL, rAN_main=0, nAN_main=0, *vAN_main=NULL;
	while (nset = bcf_sr_next_line (sr)) {
		line_main =  bcf_sr_get_line(sr, 0);
		variant_mask.push_back(false);

		//Not in the intersect of ref and main panels
		if (panels[1] && nset == 1) { n_variants_noverlap++; continue; }

		//Not a bi-alleleic variant
		if (line_main->n_allele != 2) { n_variants_multi++; continue; }

		//Unpack information if filtering
		if (filter_snp_only || (filter_min_maf > 0)) bcf_unpack(line_main, BCF_UN_STR);

		//Keep SNPs only
		if (filter_snp_only) {
			string ref = string(line_main->d.allele[0]);
			string alt = string(line_main->d.allele[1]);
			bool bref = (ref == "A") || (ref == "T") || (ref == "G") || (ref == "C");
			bool balt = (alt == "A") || (alt == "T") || (alt == "G") || (alt == "C");
			n_variants_notsnp += (!bref || !balt);
			if (!bref || !balt) continue;
		}

		//Keep common only
		if (filter_min_maf > 0) {
			rAN_main = bcf_get_info_int32(sr->readers[0].header, line_main, "AN", &vAN_main, &nAN_main);
			rAC_main = bcf_get_info_int32(sr->readers[0].header, line_main, "AC", &vAC_main, &nAC_main);
			if (nAC_main!=1) vrb.error("AC field is needed in main file for MAF filtering");
			if (nAN_main!=1) vrb.error("AN field is needed in main file for MAF filtering");
			float maf = min(vAC_main[0] * 1.0f / vAN_main[0], (vAN_main[0] - vAC_main[0]) * 1.0f / vAN_main[0]);
			n_variants_rare += (maf < filter_min_maf);
			if (maf < filter_min_maf) continue;
		}

		//Push variant information
		bcf_unpack(line_main, BCF_UN_STR);
		string chr = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
		int pos = line_main->pos + 1;
		string id = string(line_main->d.id);
		string ref = string(line_main->d.allele[0]);
		string alt = string(line_main->d.allele[1]);
		V.push(new variant (chr, pos, id, ref, alt, V.size()));


		//Flag it!
		variant_mask.back() = true;
		n_variants++;
	}
	bcf_sr_destroy(sr);
	vrb.bullet("VCF/BCF scanning done (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("#target=" + stb.str(n_main_samples) + " / #reference=" + stb.str(n_ref_samples) + " / #sites=" + stb.str(n_variants) + " / region=" + region);
	if (n_variants_noverlap) vrb.bullet(stb.str(n_variants_noverlap) + " sites removed [not in reference panel]");
	if (n_variants_multi) vrb.bullet(stb.str(n_variants_multi) + " sites removed [multi-allelic]");
	if (n_variants_notsnp) vrb.bullet(stb.str(n_variants_notsnp) + " sites removed [not SNPs]");
	if (n_variants_rare) vrb.bullet(stb.str(n_variants_rare) + " sites removed [below MAF threshold]");
	if (n_variants == 0) vrb.error("No variants to be phased!");
}

