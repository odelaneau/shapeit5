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

using std::vector;
using std::min;

void genotype_reader::scanGenotypesPlain() {
	tac.clock();
	vrb.wait("  * plain VCF/BCF scanning");

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

	//Sample processing // Needs to be improved to handle cases where sample lists do not properly overlap (in number and ordering)
	n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	int n_samples2 = bcf_hdr_nsamples(sr->readers[1].header);
	assert(n_samples == n_samples2);

	bcf1_t * line_phased = NULL, * line_unphased = NULL;
	int nset, rAC=0, nAC=0, *vAC=NULL, rAN=0, nAN=0, *vAN=NULL;
	while ((nset = bcf_sr_next_line (sr))) {
		line_unphased =  bcf_sr_get_line(sr, 0);
		line_phased =  bcf_sr_get_line(sr, 1);

		//assert(line_unphased);

		if (line_phased && line_phased->n_allele != 2) continue;
		if (line_unphased && line_unphased->n_allele != 2) continue;

		if (line_phased) {
			bcf_unpack(line_phased, BCF_UN_STR);
			std::string chr = bcf_hdr_id2name(sr->readers[1].header, line_phased->rid);
			int pos = line_phased->pos + 1;
			std::string id = std::string(line_phased->d.id);
			std::string ref = std::string(line_phased->d.allele[0]);
			std::string alt = std::string(line_phased->d.allele[1]);
			V.push(new variant (chr, pos, id, ref, alt, true, VARTYPE_SCAF));
			n_scaffold_variants++;
			n_total_variants ++;
		} else {
			bcf_unpack(line_unphased, BCF_UN_STR);
			int pos = line_unphased->pos + 1;
			if (pos >= input_start && pos <= input_stop) {
				std::string chr = bcf_hdr_id2name(sr->readers[0].header, line_unphased->rid);
				std::string id = std::string(line_unphased->d.id);
				std::string ref = std::string(line_unphased->d.allele[0]);
				std::string alt = std::string(line_unphased->d.allele[1]);
				rAN = bcf_get_info_int32(sr->readers[0].header, line_unphased, "AN", &vAN, &nAN);
				rAC = bcf_get_info_int32(sr->readers[0].header, line_unphased, "AC", &vAC, &nAC);
				assert(nAC==1 && nAN ==1);
				float maf = min(vAC[0] * 1.0f / vAN[0], (vAN[0] - vAC[0]) * 1.0f / vAN[0]);
				V.push(new variant (chr, pos, id, ref, alt, vAC[0] < (vAN[0]-vAC[0]), VARTYPE_RARE));
				n_rare_variants ++;
				n_total_variants ++;
			}
		}
	}

	bcf_sr_destroy(sr);
	if (n_rare_variants == 0) vrb.error("No variants to be phased!");
	vrb.bullet("plain VCF/BCF scanning (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_reader::scanGenotypesSparse() {
	tac.clock();
	vrb.wait("  * sparse VCF/BCF scanning");

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

	//Sample extraction
	n_samples = bcf_hdr_nsamples(sr->readers[0].header);

	//Scan files
	bcf1_t * line_phased = NULL, * line_unphased = NULL;
	int nset, rAC=0, nAC=0, *vAC=NULL, rAN=0, nAN=0, *vAN=NULL;
	while ((nset = bcf_sr_next_line (sr))) {
		line_phased =  bcf_sr_get_line(sr, 0);
		line_unphased =  bcf_sr_get_line(sr, 1);

		//Skip non-biallelic
		if (line_phased && line_phased->n_allele != 2) continue;
		if (line_unphased && line_unphased->n_allele != 2) continue;

		if (line_phased && line_unphased) vrb.error("Same record found inboth rare and common variant files");

		if (line_phased) {
			bcf_unpack(line_phased, BCF_UN_STR);
			std::string chr = bcf_hdr_id2name(sr->readers[0].header, line_phased->rid);
			int pos = line_phased->pos + 1;
			std::string id = std::string(line_phased->d.id);
			std::string ref = std::string(line_phased->d.allele[0]);
			std::string alt = std::string(line_phased->d.allele[1]);
			V.push(new variant (chr, pos, id, ref, alt, true, VARTYPE_SCAF));
			n_scaffold_variants++;
			n_total_variants ++;
		} else {
			bcf_unpack(line_unphased, BCF_UN_STR);
			int pos = line_unphased->pos + 1;
			if (pos >= input_start && pos <= input_stop) {
				std::string chr = bcf_hdr_id2name(sr->readers[1].header, line_unphased->rid);
				std::string id = std::string(line_unphased->d.id);
				std::string ref = std::string(line_unphased->d.allele[0]);
				std::string alt = std::string(line_unphased->d.allele[1]);
				rAN = bcf_get_info_int32(sr->readers[1].header, line_unphased, "AN", &vAN, &nAN);
				rAC = bcf_get_info_int32(sr->readers[1].header, line_unphased, "AC", &vAC, &nAC);
				assert(nAC==1 && nAN ==1);
				V.push(new variant (chr, pos, id, ref, alt, vAC[0] < (vAN[0]-vAC[0]), VARTYPE_RARE));
				n_rare_variants ++;
				n_total_variants ++;
			}
		}
	}

	bcf_sr_destroy(sr);
	if (n_rare_variants == 0) vrb.error("No variants to be phased!");
	vrb.bullet("VCF/BCF scanning (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

