/*******************************************************************************
 * Copyright (C) 2020 Olivier Delaneau, University of Lausanne
 * Copyright (C) 2020 Simone Rubinacci, University of Lausanne
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

#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	n_var = 0;
	n_hap = 0;
}

haplotype_set::~haplotype_set() {
	n_var = 0;
	n_hap = 0;
}

void haplotype_set::computeMAF() {
	vrb.title("Compute MAF");
	MAF = std::vector < float > (n_var, 0.0f);
	for (int v = 0 ; v < n_var ; v ++) {
		for (int h = 0 ; h < n_hap ; h++) MAF[v] += HAP[v][h];
		MAF[v] /= n_hap;
		if (MAF[v] > 0.5) MAF[v] = 1.0f - MAF[v];
	}
	vrb.bullet("Done");
}

void haplotype_set::readPositions(std::string fpos) {
	vrb.title("Reading positions to retain in [" + fpos + "]");
	std::string buffer;
	std::vector < std::string > tokens;
	input_file fd (fpos);
	while (getline(fd, buffer)) {
		stb.split(buffer, tokens);
		assert(tokens.size() == 3);
		filter_positions.insert(atoi(tokens[1].c_str()));
	}
	vrb.bullet("N=" + stb.str(filter_positions.size()));
}

bool haplotype_set::checkPos(int _pos) {
	if (filter_positions.size() == 0) return true;
	return filter_positions.count(_pos);
}

void haplotype_set::readHaplotypes(std::string fdata, std::string region, bool haploid) {
	tac.clock();

	vrb.title("Reading haplotypes in [" + fdata + "]");
	if (haploid) vrb.bullet ("Haploid mode [e.g. male chrX data]");

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, fdata.c_str());
	n_hap = (2-haploid) * bcf_hdr_nsamples(sr->readers[0].header);
	vrb.bullet ("#haps=" + stb.str(n_hap));

	for (int i = 0 ; i < n_hap ; i ++) SID.push_back(std::string(sr->readers[0].header->samples[i/(2-haploid)]));

	if (n_hap % 2) {
		n_hap --;
		SID.pop_back();
	}

	unsigned int i_variant = 0, nset = 0;
	int ngt, ngt_arr = 0, *gt_arr = NULL;
	bcf1_t * line;
	while ((nset = bcf_sr_next_line (sr))) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {

			bcf_unpack(line, BCF_UN_STR);
			std::string chr = bcf_hdr_id2name(sr->readers[0].header, line->rid);
			int pos = line->pos + 1;
			std::string id = std::string(line->d.id);
			std::string ref = std::string(line->d.allele[0]);
			std::string alt = std::string(line->d.allele[1]);

			if (checkPos(pos)) {

				CHR.push_back(chr);
				POS.push_back(pos);
				REF.push_back(ref);
				ALT.push_back(alt);
				VID.push_back(id);
				HAP.push_back(std::vector < bool > (n_hap, false));

				ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);

				//assert(ngt == n_samples);

				for(int i = 0 ; i < n_hap ; i ++) HAP.back()[i] = (bcf_gt_allele(gt_arr[i])==1);

				n_var ++;
			}
		}
	}

	free(gt_arr);
	bcf_sr_destroy(sr);

	// Report
	vrb.bullet("VCF/BCF parsing done [L="+stb.str(n_var) + " / N=" + stb.str(n_hap) + " / T=" + stb.str(tac.rel_time()*1.0/1000, 2) + "s]");
}
