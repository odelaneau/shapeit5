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

#include <io/haplotype_reader.h>

using namespace std;

haplotype_reader::haplotype_reader(haplotype_set & _H, string _region, double _minPP,  int _nthreads) : H(_H) {
	nthreads = _nthreads;
	region = _region;
	minPP = _minPP;
}

haplotype_reader::~haplotype_reader() {
	region = "";
}

void haplotype_reader::readHaplotypes(string ftruth, string festi, string ffreq, bool dupid) {
	tac.clock();
	vrb.title("Reading VCF/BCF input files");
	vrb.bullet("Validation ["  + ftruth + "]");
	vrb.bullet("Estimation ["  + festi + "]");
	vrb.bullet("Frequency  ["  + ffreq + "]");

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, ftruth.c_str()))) vrb.error("Problem opening index file for [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, festi.c_str()))) vrb.error("Problem opening index file for [" + festi + "]");
	if(!(bcf_sr_add_reader (sr, ffreq.c_str()))) vrb.error("Problem opening index file for [" + ffreq + "]");

	//IDs in truth
	int n_samples_truth = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples_truth ; i ++) {
		string sample_id = string(sr->readers[0].header->samples[i]);
		if (dupid) sample_id = sample_id + "_" + sample_id;
		H.push(sample_id);
	}
	vrb.bullet("#Validation samples = " + stb.str(n_samples_truth));

	//IDs in estimation
	int n_samples_estimated = bcf_hdr_nsamples(sr->readers[1].header);
	vector < int > mapping = vector < int > (n_samples_estimated, -1);
	for (int i = 0 ; i < n_samples_estimated ; i ++) {
		map < string, int > :: iterator itM = H.mapSamples.find(string(sr->readers[1].header->samples[i]));
		if (itM != H.mapSamples.end()) {
			mapping[i] = itM->second;
			H.IDXesti.push_back(itM->second);
		}
	}
	sort(H.IDXesti.begin(), H.IDXesti.end());
	vrb.bullet("#Estimation samples = " + stb.str(n_samples_estimated));
	vrb.bullet("#Overlapping samples = " + stb.str(H.IDXesti.size()));

	//Read GT data
	int n_variant_tot = 0;
	float * vPP = NULL;
	int rAC=0, nAC=0, *vAC = NULL, rAN=0, nAN=0, *vAN = NULL, nPP = 0, rPP = 0;
	int nset = 0, *gt_arr_t = NULL, *gt_arr_e = NULL, ngt_arr_t = 0, ngt_arr_e = 0;
	bcf1_t * line_t, * line_e, * line_f;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 3) {
			line_t =  bcf_sr_get_line(sr, 0);
			line_e =  bcf_sr_get_line(sr, 1);
			line_f =  bcf_sr_get_line(sr, 2);
			if (line_f->n_allele == 2) {
				//1. Unpack variant infos
				bcf_unpack(line_f, BCF_UN_ALL);
				H.Positions.push_back(line_f->pos + 1);
				H.RSIDs.push_back(string(line_f->d.id));
				H.REFs.push_back(string(line_f->d.allele[0]));
				H.ALTs.push_back(string(line_f->d.allele[1]));

				rAN = bcf_get_info_int32(sr->readers[2].header, line_f, "AN", &vAN, &nAN);
				rAC = bcf_get_info_int32(sr->readers[2].header, line_f, "AC", &vAC, &nAC);
				if (nAC!=1) vrb.error("AC field is needed");
				if (nAN!=1) vrb.error("AN field is needed");
				H.MAC.push_back(min(vAC[0], (vAN[0] - vAC[0])));
				H.MinorAlleles.push_back(vAC[0] < (vAN[0] - vAC[0]));

				//2. Validation
				bcf_get_genotypes(sr->readers[0].header, line_t, &gt_arr_t, &ngt_arr_t);
				for(int h = 0 ; h < 2 * n_samples_truth ; h += 2) {
					bool a0 = bcf_gt_allele(gt_arr_t[h+0])!=0;
					bool a1 = bcf_gt_allele(gt_arr_t[h+1])!=0;
					bool mi = (gt_arr_t[h+0] == bcf_gt_missing || gt_arr_t[h+1] == bcf_gt_missing);
					H.Htrue[h+0].push_back(a0);
					H.Htrue[h+1].push_back(a1);
					H.Missing[h/2].push_back(gt_arr_t[h+0] == bcf_gt_missing || gt_arr_t[h+1] == bcf_gt_missing);
				}

				//3. Estimation
				bcf_get_genotypes(sr->readers[1].header, line_e, &gt_arr_e, &ngt_arr_e);
				for(int h = 0 ; h < 2 * n_samples_estimated ; h += 2) {
					int index = mapping[h/2];
					if (index >= 0) {
						bool a0 = bcf_gt_allele(gt_arr_e[h+0])!=0;
						bool a1 = bcf_gt_allele(gt_arr_e[h+1])!=0;
						H.Hesti[2*index+0].push_back(a0);
						H.Hesti[2*index+1].push_back(a1);
					}
				}

				//4. Probabilities
				rPP = bcf_get_format_float(sr->readers[1].header, line_e, "PP", &vPP, &nPP);
				if (rPP == n_samples_estimated) {
					for(int i = 0 ; i < n_samples_estimated ; i ++) {
						int index = mapping[i];
						if (index >= 0) {
							H.Hprob[index].push_back(!bcf_float_is_missing(vPP[i]));
							H.Estimated[index].push_back(true);
							if (!bcf_float_is_missing(vPP[i])) {
								string key = stb.str(H.Hprob[index].size() - 1) + "_" + stb.str(index);
								H.Vprob.insert(pair < string, float > ( key, vPP[i]));
								if (vPP[i] <= minPP) H.Estimated[index].back() = false;
							}
						}
					}
				} else {
					for(int i = 0 ; i < n_samples_estimated ; i ++) {
						int index = mapping[i];
						if (index >= 0) {
							H.Hprob[index].push_back(false);
							H.Estimated[index].push_back(true);
						}
					}
				}

				H.n_variants ++;

			}
		}
		n_variant_tot ++;
		if (n_variant_tot % 10000 == 0) vrb.bullet (stb.str(n_variant_tot) + " lines processed");
	}
	vrb.bullet("#Total variants = " + stb.str(n_variant_tot));
	vrb.bullet("#Overlapping variants = " + stb.str(H.n_variants));
	vrb.bullet("#Prob stored [PP field] = " + stb.str(H.Vprob.size()));
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
	free(gt_arr_t); free(gt_arr_e);
	bcf_sr_destroy(sr);
}


void haplotype_reader::readHaplotypes(string ftruth, string festi, bool dupid) {
	tac.clock();
	vrb.title("Reading VCF/BCF input files");
	vrb.bullet("Validation ["  + ftruth + "]");
	vrb.bullet("Estimation ["  + festi + "]");

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, ftruth.c_str()))) vrb.error("Problem opening index file for [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, festi.c_str()))) vrb.error("Problem opening index file for [" + festi + "]");

	//IDs in truth
	int n_samples_truth = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples_truth ; i ++) {
		string sample_id = string(sr->readers[0].header->samples[i]);
		if (dupid) sample_id = sample_id + "_" + sample_id;
		H.push(sample_id);
	}
	vrb.bullet("#Validation samples = " + stb.str(n_samples_truth));

	//IDs in estimation
	int n_samples_estimated = bcf_hdr_nsamples(sr->readers[1].header);
	vector < int > mapping = vector < int > (n_samples_estimated, -1);
	for (int i = 0 ; i < n_samples_estimated ; i ++) {
		map < string, int > :: iterator itM = H.mapSamples.find(string(sr->readers[1].header->samples[i]));
		if (itM != H.mapSamples.end()) {
			mapping[i] = itM->second;
			H.IDXesti.push_back(itM->second);
		}
	}
	sort(H.IDXesti.begin(), H.IDXesti.end());
	vrb.bullet("#Estimation samples = " + stb.str(n_samples_estimated));
	vrb.bullet("#Overlapping samples = " + stb.str(H.IDXesti.size()));

	//Read GT data
	int n_variant_tot = 0;
	float * vPP = NULL;
	int rAC=0, nAC=0, *vAC = NULL, rAN=0, nAN=0, *vAN = NULL, nPP = 0, rPP = 0;
	int nset = 0, *gt_arr_t = NULL, *gt_arr_e = NULL, ngt_arr_t = 0, ngt_arr_e = 0;
	bcf1_t * line_t, * line_e;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_t =  bcf_sr_get_line(sr, 0);
			line_e =  bcf_sr_get_line(sr, 1);
			if (line_t->n_allele == 2) {
				//1. Unpack variant infos
				bcf_unpack(line_t, BCF_UN_ALL);
				H.Positions.push_back(line_t->pos + 1);
				H.RSIDs.push_back(string(line_t->d.id));
				H.REFs.push_back(string(line_t->d.allele[0]));
				H.ALTs.push_back(string(line_t->d.allele[1]));

				//2. Validation
				int vAC = 0, vAN = 0;
				bcf_get_genotypes(sr->readers[0].header, line_t, &gt_arr_t, &ngt_arr_t);
				for(int h = 0 ; h < 2 * n_samples_truth ; h += 2) {
					bool a0 = bcf_gt_allele(gt_arr_t[h+0])!=0;
					bool a1 = bcf_gt_allele(gt_arr_t[h+1])!=0;
					bool mi = (gt_arr_t[h+0] == bcf_gt_missing || gt_arr_t[h+1] == bcf_gt_missing);
					H.Htrue[h+0].push_back(a0);
					H.Htrue[h+1].push_back(a1);
					H.Missing[h/2].push_back(gt_arr_t[h+0] == bcf_gt_missing || gt_arr_t[h+1] == bcf_gt_missing);
					if (!mi) {
						vAC += a0+a1;
						vAN += 2;
					}
				}
				H.MAC.push_back(min(vAC, (vAN - vAC)));
				H.MinorAlleles.push_back(vAC < (vAN - vAC));

				//3. Estimation
				bcf_get_genotypes(sr->readers[1].header, line_e, &gt_arr_e, &ngt_arr_e);
				for(int h = 0 ; h < 2 * n_samples_estimated ; h += 2) {
					int index = mapping[h/2];
					if (index >= 0) {
						bool a0 = bcf_gt_allele(gt_arr_e[h+0])!=0;
						bool a1 = bcf_gt_allele(gt_arr_e[h+1])!=0;
						H.Hesti[2*index+0].push_back(a0);
						H.Hesti[2*index+1].push_back(a1);
					}
				}

				//4. Probabilities
				rPP = bcf_get_format_float(sr->readers[1].header, line_e, "PP", &vPP, &nPP);
				if (rPP == n_samples_estimated) {
					for(int i = 0 ; i < n_samples_estimated ; i ++) {
						int index = mapping[i];
						if (index >= 0) {
							H.Hprob[index].push_back(!bcf_float_is_missing(vPP[i]));
							H.Estimated[index].push_back(true);
							if (!bcf_float_is_missing(vPP[i])) {
								string key = stb.str(H.Hprob[index].size() - 1) + "_" + stb.str(index);
								H.Vprob.insert(pair < string, float > ( key, vPP[i]));
								if (vPP[i] <= minPP) H.Estimated[index].back() = false;
							}
						}
					}
				} else {
					for(int i = 0 ; i < n_samples_estimated ; i ++) {
						int index = mapping[i];
						if (index >= 0) {
							H.Hprob[index].push_back(false);
							H.Estimated[index].push_back(true);
						}
					}
				}

				H.n_variants ++;
			}
		}
		n_variant_tot ++;
		if (n_variant_tot % 10000 == 0) vrb.bullet (stb.str(n_variant_tot) + " lines processed");
	}
	vrb.bullet("#Total variants = " + stb.str(n_variant_tot));
	vrb.bullet("#Overlapping variants = " + stb.str(H.n_variants));
	vrb.bullet("#Prob stored [PP field] = " + stb.str(H.Vprob.size()));
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
	free(gt_arr_t); free(gt_arr_e);
	bcf_sr_destroy(sr);
}
