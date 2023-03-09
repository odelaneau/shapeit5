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

#include "../../versions/versions.h"

#include <io/sparse2plain.h>
#include <objects/rare_genotype.h>

using namespace std;

sparse2plain::sparse2plain(string _sparse_common, string _sparse_rare, string _plain_bcf, string _region, int _nthreads) {
	file_full_bcf = _plain_bcf;
	file_comm_bcf = _sparse_common;
	file_rare_bcf = _sparse_rare;
	file_rare_bin = stb.get_name_from_vcf(_sparse_rare) + ".bin";
	file_rare_prb = stb.get_name_from_vcf(_sparse_rare) + ".prb";
	region = _region;
	nthreads = _nthreads;

	vector < string > tokens;
	stb.split(region, tokens, ":");
	if (tokens.size() == 1 || tokens.size() == 2) contig = tokens[0];
	else vrb.error("Could not parse --region");
}

sparse2plain::~sparse2plain() {
}

void sparse2plain::convert() {
	tac.clock();
	vrb.title("Converting from sparse BCF file to plain BCF");
	vrb.bullet("Region        : " + region);
	vrb.bullet("Contig        : " + contig);

	//Opening sparse VCF/BCF file
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "]");
	if (!(bcf_sr_add_reader (sr, file_comm_bcf.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + file_comm_bcf + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + file_comm_bcf + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + file_comm_bcf + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + file_comm_bcf + "]: unknown error"); break;
		}
	}
	if (!(bcf_sr_add_reader (sr, file_rare_bcf.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + file_rare_bcf + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + file_rare_bcf + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + file_rare_bcf + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + file_rare_bcf + "]: unknown error"); break;
		}
	}

	//Checking binary files to open
	int id = bcf_hdr_id2int(sr->readers[1].header, BCF_DT_ID, "SGEN");
	if (!bcf_hdr_idinfo_exists(sr->readers[1].header, BCF_HL_INFO, id)) vrb.error("Cannot retrieve SGEN flag in VCF/BCF for sparse rare variants.");
	id = bcf_hdr_id2int(sr->readers[1].header, BCF_DT_ID, "SPRB");
	bool hasProbs = bcf_hdr_idinfo_exists(sr->readers[1].header, BCF_HL_INFO, id);
	if (hasProbs) vrb.bullet("Sparse probability file detected");

	//Opening sparse BIN file
	std::ifstream fp_rare_bin (file_rare_bin, std::ios::in | std::ios::binary);
	if (!fp_rare_bin) vrb.error("Cannot open " + file_rare_bin + " for reading, check permissions");
	std::ifstream fp_rare_prb;
	if (hasProbs) {
		fp_rare_prb.open(file_rare_prb, std::ios::in | std::ios::binary);
		if (!fp_rare_prb) vrb.error("Cannot open " + file_rare_prb + " for reading, check permissions");
	}

	//Opening plain file
	htsFile * fp_full_bcf = hts_open(file_full_bcf.c_str(), "wb");
	if (!fp_full_bcf) vrb.error("Cannot open " + file_full_bcf + " for writing, check permissions");
	if (nthreads > 1) hts_set_threads(fp_full_bcf, nthreads);

	//Create and write VCF header for plain version
	bcf_hdr_t * hdr_full_bcf = bcf_hdr_init("w");
	bcf_hdr_append(hdr_full_bcf, string("##source=shapeit5 convert v" + string(SCFTLS_VERSION)).c_str());
	bcf_hdr_append(hdr_full_bcf, string("##contig=<ID="+ contig + ">").c_str());
	bcf_hdr_append(hdr_full_bcf, "##INFO=<ID=AN,Number=1,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr_full_bcf, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr_full_bcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	if (hasProbs) bcf_hdr_append(hdr_full_bcf, "##FORMAT=<ID=PP,Number=1,Type=Float,Description=\"Phasing confidence\">");
	int n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples ; i ++) bcf_hdr_add_sample(hdr_full_bcf, sr->readers[0].header->samples[i]);
	bcf_hdr_add_sample(hdr_full_bcf, NULL);
	if (bcf_hdr_write(fp_full_bcf, hdr_full_bcf) < 0) vrb.error("Failing to write VCF/header for common variants");

	//
	bcf1_t * line_input_comm = NULL;
	bcf1_t * line_input_rare = NULL;
	bcf1_t * line_output = bcf_init1();
	int ngt = 0, rgt = 0, *vgt = NULL;
	int nac = 0, rac = 0, *vac = NULL;
	int nan = 0, ran = 0, *van = NULL;
	int nsk = 0, rsk = 0, *vsk = NULL;
	int ngs = 0, rgs = 0, *vgs = (int*)malloc(n_samples*2*sizeof(int));
	float * probabilities = (float*)malloc(n_samples*sizeof(float));
	unsigned int * rg_buffer = (unsigned int *)malloc(n_samples * sizeof(unsigned int));
	float * rp_buffer = (float*)malloc(n_samples*sizeof(float));

	//
	unsigned int nrare = 0, ncomm = 0, nfull = 0, nset = 0;
	while ((nset = bcf_sr_next_line (sr))) {
		line_input_comm = bcf_sr_get_line(sr, 0);
		line_input_rare = bcf_sr_get_line(sr, 1);

		//Variant information
		string chr, id, ref, alt;
		int pos;

		//Skip not bi-allelic
		if (line_input_comm && line_input_comm->n_allele != 2) continue;
		if (line_input_rare && line_input_rare->n_allele != 2) continue;

		//
		if (nset == 2) vrb.error("Duplicate variant at pos= " + stb.str(pos+1));

		//Reading in variant information for COMMON
		if (line_input_comm) {
			chr = bcf_hdr_id2name(sr->readers[0].header, line_input_comm->rid);
			pos = line_input_comm->pos;
			id = string(line_input_comm->d.id);
			ref = string(line_input_comm->d.allele[0]);
			alt = string(line_input_comm->d.allele[1]);

			ran = bcf_get_info_int32(sr->readers[0].header, line_input_comm, "AN", &van, &nan); if (nan!=1) vrb.error("AN field is needed in common file " + stb.str(nan));
			rac = bcf_get_info_int32(sr->readers[0].header, line_input_comm, "AC", &vac, &nac); if (nac!=1) vrb.error("AC field is needed in common file");

			rgt = bcf_get_genotypes(sr->readers[0].header, line_input_comm, &vgt, &ngt);
			ncomm++;
		}

		//Reading in variant information for RARE
		if (line_input_rare) {
			chr = bcf_hdr_id2name(sr->readers[1].header, line_input_rare->rid);
			pos = line_input_rare->pos;
			id = string(line_input_rare->d.id);
			ref = string(line_input_rare->d.allele[0]);
			alt = string(line_input_rare->d.allele[1]);

			//Process MAF
			ran = bcf_get_info_int32(sr->readers[1].header, line_input_rare, "AN", &van, &nan); if (nan!=1) vrb.error("AN field is needed in rare file");
			rac = bcf_get_info_int32(sr->readers[1].header, line_input_rare, "AC", &vac, &nac); if (nac!=1) vrb.error("AC field is needed in rare file");
			float curraf = vac[0] * 1.0f / van[0];
			for (int i = 0 ; i < 2*n_samples; i ++) vgs[i] = bcf_gt_phased(curraf >= 0.5f);
			if (hasProbs) for (int i = 0 ; i < n_samples; i ++) bcf_float_set_missing(probabilities[i]);
			ngs = 2*n_samples;

			//Process sparse genotypes addressing
			rsk = bcf_get_info_int32(sr->readers[1].header, line_input_rare, "SGEN", &vsk, &nsk); if (nsk!=3) vrb.error("SEEK field is needed in rare file");
			uint64_t sgenotypes = vsk[0];
			sgenotypes *= MOD30BITS;
			sgenotypes += vsk[1];
			uint32_t ngenotypes = vsk[2];

			//Read sparse genotypes in file
			fp_rare_bin.seekg(sgenotypes * sizeof(unsigned int), fp_rare_bin.beg);
			fp_rare_bin.read((char *)rg_buffer, ngenotypes * sizeof(unsigned int));

			//Read sparse probabilities
			if (hasProbs) {
				fp_rare_prb.seekg(sgenotypes * sizeof(float), fp_rare_prb.beg);
				fp_rare_prb.read((char *)rp_buffer, ngenotypes * sizeof(float));
			}

			//Parse sparse genotypes and probs
			for (int r = 0 ; r < ngenotypes ; r++) {
				//
				rare_genotype rg;
				rg.set(rg_buffer[r]);
				if (rg.mis) {
					vgs[2*rg.idx+0] = bcf_gt_missing;
					vgs[2*rg.idx+1] = bcf_gt_missing;
				} else if (rg.het) {
					if (rg.pha) {
						vgs[2*rg.idx+0] = bcf_gt_phased(rg.al0);
						vgs[2*rg.idx+1] = bcf_gt_phased(rg.al1);
					} else {
						vgs[2*rg.idx+0] = bcf_gt_unphased(rg.al0);
						vgs[2*rg.idx+1] = bcf_gt_unphased(rg.al1);
					}
				}
				//
				probabilities[rg.idx] = rp_buffer[r];
			}

			//
			nrare++;
		}

		//Writing variant information
		bcf_clear1(line_output);
		line_output->rid = bcf_hdr_name2id(hdr_full_bcf, chr.c_str());
		line_output->pos = pos;
		bcf_update_id(hdr_full_bcf, line_output, id.c_str());
		string alleles = ref + "," + alt;
		bcf_update_alleles_str(hdr_full_bcf, line_output, alleles.c_str());

		//Writing genotypes
		if (line_input_comm) bcf_update_genotypes(hdr_full_bcf, line_output, vgt, ngt);
		if (line_input_rare) {
			bcf_update_genotypes(hdr_full_bcf, line_output, vgs, ngs);
			if (hasProbs) bcf_update_format_float(hdr_full_bcf, line_output, "PP", probabilities, bcf_hdr_nsamples(hdr_full_bcf));
		}
		bcf_update_info_int32(hdr_full_bcf, line_output, "AC", vac, 1);
		bcf_update_info_int32(hdr_full_bcf, line_output, "AN", van, 1);
		if (bcf_write1(fp_full_bcf, hdr_full_bcf, line_output) < 0) vrb.error("Failing to write VCF/record");
		nfull ++;

		if (nfull % 10000 == 0) vrb.bullet("VCF/BCF writing [nfull=" + stb.str(nfull) + ", ncomm=" + stb.str(ncomm) + ", nrare=" + stb.str(nrare) + "]");
	}

	//Closing stuffs
	free(rg_buffer);
	free(rp_buffer);
	free(probabilities);
	free(vgt);
	free(vgs);
	free(vac);
	free(van);
	free(vsk);
	bcf_sr_destroy(sr);
	fp_rare_bin.close();
	fp_rare_prb.close();
	if (hts_close(fp_full_bcf)) vrb.error("Non zero status when closing VCF/BCF file descriptor");

	// Report
	vrb.bullet("VCF/BCF writing done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	vrb.bullet("Indexing ["+file_full_bcf + "]");
	if (bcf_index_build3(file_full_bcf.c_str(), NULL, 14, nthreads) < 0) vrb.error("Fail to index file");
}

