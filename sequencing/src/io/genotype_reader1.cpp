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

genotype_reader::genotype_reader(haplotype_set & _H, genotype_set & _G, variant_map & _V, string _region, bool _use_PS_field, int _nthreads) : H(_H), G(_G), V(_V) {
	nthreads = _nthreads;
	n_common_variants = 0;
	n_rare_variants = 0;
	n_target_samples = 0;
	n_reference_samples = 0;
	region = _region;
	n_genotypes = vector < unsigned long > (7,0);
	use_PS_field = _use_PS_field;
	maf_threshold = RARE_VARIANT_FREQ;
}

genotype_reader::~genotype_reader() {
	n_common_variants = 0;
	n_rare_variants = 0;
	n_target_samples = 0;
	n_reference_samples = 0;
	region = "";
}

void genotype_reader::allocateGenotypes() {
	assert(n_common_variants != 0 && (n_target_samples+n_reference_samples) != 0);
	G.allocate(n_target_samples, n_common_variants);
	H.allocate(n_target_samples, n_reference_samples, n_common_variants, n_rare_variants, flag_common);
}

void genotype_reader::setPScodes(int * ps_arr, int nps) {
	if (nps != n_target_samples) PScodes.clear();
	else {
		int ps_idx = 0;
		unordered_map < int , int > :: iterator it;
		PScodes = vector < int > (n_target_samples, 0);
		for (int i = 0 ; i < n_target_samples ; i ++) {
			if (ps_arr[i] != bcf_int32_missing) {
				it = PSmap.find(ps_arr[i]);
				if (it == PSmap.end()) {
					ps_idx = PSmap.size()+1;
					PSmap.insert(pair < int, int > (ps_arr[i], ps_idx));
				} else ps_idx = it->second;
				PScodes[i] = ps_idx;
			}
		}
	}
}

void genotype_reader::scanGenotypes(string fmain) {
	vrb.wait("  * VCF/BCF scanning");
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fmain.c_str()))) vrb.error("Problem opening index file for [" + fmain + "]");

	n_common_variants = 0;
	n_rare_variants = 0;
	n_total_variants = 0;

	n_target_samples = bcf_hdr_nsamples(sr->readers[0].header);
	n_reference_samples = 0;
	if (n_target_samples < 50) vrb.error("Statistical phasing for less than 50 individuals is not permitted, use a reference panel to solve this issue!");

	bcf1_t * line;
	int rAC=0, nAC=0, *vAC = NULL;
	int rAN=0, nAN=0, *vAN = NULL;

	while(bcf_sr_next_line (sr)) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {
			bcf_unpack(line, BCF_UN_INFO);
			string chr = bcf_hdr_id2name(sr->readers[0].header, line->rid);
			unsigned int pos = line->pos + 1;
			string id = string(line->d.id);
			string ref = string(line->d.allele[0]);
			string alt = string(line->d.allele[1]);
			variant * new_variant = new variant (chr, pos, id, ref, alt);
			rAC = bcf_get_info_int32(sr->readers[0].header, line, "AC", &vAC, &nAC);
			rAN = bcf_get_info_int32(sr->readers[0].header, line, "AN", &vAN, &nAN);
			if ((nAC!=1)||(nAN!=1)) vrb.error("AC/AN INFO fields are needed");
			new_variant->minor = (vAC[0] < (vAN[0] - vAC[0]));

			float current_maf = min(vAC[0] * 1.0 / vAN[0], (vAN[0] - vAC[0]) * 1.0 / vAN[0]);
			if (current_maf < maf_threshold) {
				flag_common.push_back(false);
				new_variant->rare = true;
				n_rare_variants++;
			} else {
				flag_common.push_back(true);
				new_variant->rare = false;
				n_common_variants++;
			}
			n_total_variants++;
			V.push(new_variant);
		}
	}
	if (flag_common.size() == 0) flag_common = vector < bool > (n_total_variants, true);
	bcf_sr_destroy(sr);
	if (n_common_variants == 0) vrb.error("No variants to be phased in [" + fmain + "]");
	if (n_rare_variants == 0)
		vrb.bullet("VCF/BCF scanning [N=" + stb.str(n_target_samples) + " / L=" + stb.str(n_common_variants) + " / Reg=" + region + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	else
		vrb.bullet("VCF/BCF scanning [N=" + stb.str(n_target_samples) + " / Lc=" + stb.str(n_common_variants)  + " / Lr=" + stb.str(n_rare_variants) + " / Reg=" + region + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_reader::scanGenotypes(string fmain, string fref) {
	vrb.wait("  * VCF/BCF scanning");
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fmain.c_str()))) vrb.error("Problem opening index file for [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fref.c_str()))) vrb.error("Problem opening index file for [" + fref + "]");

	n_common_variants = 0;
	n_rare_variants = 0;
	n_total_variants = 0;

	n_target_samples = bcf_hdr_nsamples(sr->readers[0].header);
	n_reference_samples = bcf_hdr_nsamples(sr->readers[1].header);
	if ((n_reference_samples + n_target_samples) < 50) vrb.error("Statistical phasing for less than 50 individuals is not permitted!");

	int nset;
	bcf1_t * line_target, * line_reference;
	int rAC=0, nAC=0, *vAC = NULL;
	int rAN=0, nAN=0, *vAN = NULL;

	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_target =  bcf_sr_get_line(sr, 0);
			line_reference =  bcf_sr_get_line(sr, 1);
			if (line_target->n_allele == 2 && line_reference->n_allele == 2) {
				bcf_unpack(line_reference, BCF_UN_INFO);
				bcf_unpack(line_target, BCF_UN_INFO);
				string chr = bcf_hdr_id2name(sr->readers[1].header, line_reference->rid);
				unsigned int pos = line_reference->pos + 1;
				string id = string(line_reference->d.id);
				string ref = string(line_reference->d.allele[0]);
				string alt = string(line_reference->d.allele[1]);
				variant * new_variant = new variant (chr, pos, id, ref, alt);
				rAC = bcf_get_info_int32(sr->readers[1].header, line_reference, "AC", &vAC, &nAC);
				rAN = bcf_get_info_int32(sr->readers[1].header, line_reference, "AN", &vAN, &nAN);
				if ((nAC!=1)||(nAN!=1)) vrb.error("AC/AN INFO fields are needed");
				new_variant->minor = (vAC[0] < (vAN[0] - vAC[0]));
				float current_maf = min(vAC[0] * 1.0 / vAN[0], (vAN[0] - vAC[0]) * 1.0 / vAN[0]);
				if (current_maf < maf_threshold) {
					flag_common.push_back(false);
					new_variant->rare = true;
					n_rare_variants++;
				} else {
					flag_common.push_back(true);
					new_variant->rare = false;
					n_common_variants++;
				}
				n_total_variants++;
				V.push(new_variant);
			}
		}
	}

	bcf_sr_destroy(sr);
	if (n_common_variants == 0) vrb.error("No variants to be phased in files");
	vrb.bullet("VCF/BCF scanning [Nm=" + stb.str(n_target_samples) + " / Nr=" + stb.str(n_reference_samples) + " / L=" + stb.str(n_common_variants) + " / Reg=" + region + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
