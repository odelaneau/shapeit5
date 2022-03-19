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
	if (!(bcf_sr_add_reader (sr, funphased.c_str()))) vrb.error("Problem opening index file for [" + funphased + "]");
	if (!(bcf_sr_add_reader (sr, fphased.c_str()))) vrb.error("Problem opening index file for [" + fphased + "]");

	//Sample processing // Needs to be improved to handle cases where sample lists do not properly overlap (in number and ordering)
	n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	int n_samples2 = bcf_hdr_nsamples(sr->readers[1].header);
	assert(n_samples == n_samples2);

	bcf1_t * line_phased, * line_unphased;
	int nset, rAC=0, nAC=0, *vAC=NULL, rAN=0, nAN=0, *vAN=NULL;
	while (nset = bcf_sr_next_line (sr)) {
		line_unphased =  bcf_sr_get_line(sr, 0);
		line_phased =  bcf_sr_get_line(sr, 1);

		if (line_phased && line_phased->n_allele != 2) continue;
		if (line_unphased && line_phased->n_allele != 2) continue;

		if (line_phased) {
			bcf_unpack(line_phased, BCF_UN_STR);
			string chr = bcf_hdr_id2name(sr->readers[1].header, line_phased->rid);
			int pos = line_phased->pos + 1;
			string id = string(line_phased->d.id);
			string ref = string(line_phased->d.allele[0]);
			string alt = string(line_phased->d.allele[1]);
			V.push(new variant (chr, pos, id, ref, alt, true, VARTYPE_SCAF));
			n_scaffold_variants++;
		} else {
			bcf_unpack(line_unphased, BCF_UN_STR);
			string chr = bcf_hdr_id2name(sr->readers[0].header, line_unphased->rid);
			int pos = line_unphased->pos + 1;
			string id = string(line_unphased->d.id);
			string ref = string(line_unphased->d.allele[0]);
			string alt = string(line_unphased->d.allele[1]);
			rAN = bcf_get_info_int32(sr->readers[0].header, line_unphased, "AN", &vAN, &nAN);
			rAC = bcf_get_info_int32(sr->readers[0].header, line_unphased, "AC", &vAC, &nAC);
			assert(nAC==1 && nAN ==1);
			float maf = min(vAC[0] * 1.0f / vAN[0], (vAN[0] - vAC[0]) * 1.0f / vAN[0]);
			V.push(new variant (chr, pos, id, ref, alt, vAC[0] < (vAN[0]-vAC[0]), (maf<minmaf)?VARTYPE_RARE:VARTYPE_COMM));
			n_rare_variants += (maf<minmaf);
			n_common_variants += (maf>=minmaf);
		}

		n_total_variants ++;
	}

	bcf_sr_destroy(sr);
	if ((n_rare_variants + n_common_variants) == 0) vrb.error("No variants to be phased!");
	vrb.bullet("VCF/BCF scanning done (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

