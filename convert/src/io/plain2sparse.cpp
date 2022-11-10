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

#include <io/plain2sparse.h>
#include <objects/rare_genotype.h>

using namespace std;

plain2sparse::plain2sparse(string _plain_vcf, string _sparse_prefix, string _region, int _nthreads, float _minmaf) {
	file_full_vcf = _plain_vcf;
	file_comm_bcf = _sparse_prefix +".comm.bcf";
	file_rare_bcf = _sparse_prefix +".rare.bcf";
	file_rare_bin = _sparse_prefix +".rare.bin";
	region = _region;
	minmaf = _minmaf;
	nthreads = _nthreads;

	vector < string > tokens;
	stb.split(region, tokens, ":");
	if (tokens.size() == 1 || tokens.size() == 2) contig = tokens[0];
	else vrb.error("Could not parse --region");
}

plain2sparse::~plain2sparse() {
}

void plain2sparse::convert() {
	tac.clock();
	vrb.title("Converting from plain VCF/BCF file to sparse BCF");
	vrb.bullet("Region        : " + region);
	vrb.bullet("Contig        : " + contig);
	vrb.bullet("MAF threshold : " + stb.str(minmaf));

	//Opening plain VCF/BCF file
	bcf_srs_t * sr =  bcf_sr_init();
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "]");
	if (!(bcf_sr_add_reader (sr, file_full_vcf.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + file_full_vcf + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + file_full_vcf + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + file_full_vcf + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + file_full_vcf + "]: unknown error"); break;
		}
	}

	//Opening all sparse files
	htsFile * fp_comm_bcf = hts_open(file_comm_bcf.c_str(), "wb");
	if (!fp_comm_bcf) vrb.error("Cannot open " + file_comm_bcf + " for writing, check permissions");
	if (nthreads > 1) hts_set_threads(fp_comm_bcf, nthreads);
	htsFile * fp_rare_bcf = hts_open(file_rare_bcf.c_str(), "wb");
	if (!fp_rare_bcf) vrb.error("Cannot open " + file_rare_bcf + " for writing, check permissions");
	ofstream fp_rare_bin (file_rare_bin, std::ios::out | std::ios::binary);
	if (!fp_rare_bin) vrb.error("Cannot open " + file_rare_bin + " for writing, check permissions");

	//Create and write VCF header for common
	bcf_hdr_t * hdr_comm_bcf = bcf_hdr_init("w");
	bcf_hdr_append(hdr_comm_bcf, string("##source=shapeit5 convert v" + string(CONVER_VERSION)).c_str());
	bcf_hdr_append(hdr_comm_bcf, string("##contig=<ID="+ contig + ">").c_str());
	bcf_hdr_append(hdr_comm_bcf, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr_comm_bcf, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr_comm_bcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	int n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples ; i ++) bcf_hdr_add_sample(hdr_comm_bcf, sr->readers[0].header->samples[i]);
	bcf_hdr_add_sample(hdr_comm_bcf, NULL);
	if (bcf_hdr_write(fp_comm_bcf, hdr_comm_bcf) < 0) vrb.error("Failing to write VCF/header for common variants");

	//Create and write VCF header for rare
	bcf_hdr_t * hdr_rare_bcf = bcf_hdr_init("w");
	bcf_hdr_append(hdr_rare_bcf, string("##source=shapeit5 convert v" + string(CONVER_VERSION)).c_str());
	bcf_hdr_append(hdr_rare_bcf, string("##contig=<ID="+ contig + ">").c_str());
	bcf_hdr_append(hdr_rare_bcf, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr_rare_bcf, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr_rare_bcf, "##INFO=<ID=SEEK,Number=2,Type=Integer,Description=\"Index of the first record and number of records for the variant site in associated bin file\">");
	bcf_hdr_add_sample(hdr_rare_bcf, NULL);
	if (bcf_hdr_write(fp_rare_bcf, hdr_rare_bcf) < 0) vrb.error("Failing to write VCF/header for rare variants");


	//
	bcf1_t * line_input = NULL;
	bcf1_t * line_output_comm = bcf_init1();
	bcf1_t * line_output_rare = bcf_init1();
	int ngt = 0, rgt = 0, *vgt = NULL;
	int nac = 0, rac = 0, *vac = NULL;
	int nan = 0, ran = 0, *van = NULL;
	int nsk, rsk, *vsk = (int*)malloc(2*sizeof(int));

	//
	unsigned int nseek = 0, nrare = 0, ncomm = 0, nfull = 0;
	while (bcf_sr_next_line (sr)) {
		line_input = bcf_sr_get_line(sr, 0);

		//Skip not bi-allelic
		if (line_input->n_allele != 2) continue;

		//Get variant infos
		string chr = bcf_hdr_id2name(sr->readers[0].header, line_input->rid);
		int pos = line_input->pos;
		string id = string(line_input->d.id);
		string ref = string(line_input->d.allele[0]);
		string alt = string(line_input->d.allele[1]);

		//Check variant MAF
		ran = bcf_get_info_int32(sr->readers[0].header, line_input, "AN", &van, &nan);
		rac = bcf_get_info_int32(sr->readers[0].header, line_input, "AC", &vac, &nac);
		if (nan!=1) vrb.error("AN field is needed in main file for MAF filtering");
		if (nac!=1) vrb.error("AC field is needed in main file for MAF filtering");
		float currmaf = min(vac[0] * 1.0f / van[0], (van[0] - vac[0]) * 1.0f / van[0]);

		//Get genotypes
		rgt = bcf_get_genotypes(sr->readers[0].header, line_input, &vgt, &ngt);

		//RARE VARIANT
		if (currmaf < minmaf) {
			bcf_clear1(line_output_rare);
			line_output_rare->rid = bcf_hdr_name2id(hdr_rare_bcf, chr.c_str());
			line_output_rare->pos = pos;
			bcf_update_id(hdr_rare_bcf, line_output_rare, id.c_str());
			string alleles = ref + "," + alt;
			bcf_update_alleles_str(hdr_rare_bcf, line_output_rare, alleles.c_str());

			//

			vsk[0] = nseek; vsk[1] = 0;
			bool minor_allele = ((van[0] - vac[0]) > vac[0]);
			for(int i = 0 ; i < 2 * n_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(vgt[i+0])==1);
				bool a1 = (bcf_gt_allele(vgt[i+1])==1);
				bool mi = (vgt[i+0] == bcf_gt_missing || vgt[i+1] == bcf_gt_missing);
				bool ph = (bcf_gt_is_phased(vgt[i+0]) || bcf_gt_is_phased(vgt[i+1]));
				if ( a0 == minor_allele || a1 == minor_allele || mi) {
					rare_genotype rg_struct = rare_genotype(i/2, (a0!=a1), mi, a0, a1, ph);
					unsigned int rg_int = rg_struct.get();
					fp_rare_bin.write(reinterpret_cast < char * > (&rg_int), sizeof(unsigned int));
					vsk[1] ++;
					nseek ++;
				}
			}

			//
			bcf_update_info_int32(hdr_rare_bcf, line_output_rare, "SEEK", vsk, 2);
			bcf_update_info_int32(hdr_rare_bcf, line_output_rare, "AC", vac, 1);
			bcf_update_info_int32(hdr_rare_bcf, line_output_rare, "AN", van, 1);
			if (bcf_write1(fp_rare_bcf, hdr_rare_bcf, line_output_rare) < 0) vrb.error("Failing to write VCF/record for rare variants");
			nrare++;
		}

		//COMMON VARIANT
		if (currmaf >= minmaf) {
			bcf_clear1(line_output_comm);
			line_output_comm->rid = bcf_hdr_name2id(hdr_comm_bcf, chr.c_str());
			line_output_comm->pos = pos;
			bcf_update_id(hdr_comm_bcf, line_output_comm, id.c_str());
			string alleles = ref + "," + alt;
			bcf_update_alleles_str(hdr_comm_bcf, line_output_comm, alleles.c_str());
			bcf_update_info_int32(hdr_comm_bcf, line_output_comm, "AC", vac, 1);
			bcf_update_info_int32(hdr_comm_bcf, line_output_comm, "AN", van, 1);
			bcf_update_genotypes(hdr_comm_bcf, line_output_comm, vgt, ngt);
			if (bcf_write1(fp_comm_bcf, hdr_comm_bcf, line_output_comm) < 0) vrb.error("Failing to write VCF/record for common variants");
			ncomm ++;
		}
		nfull ++;

		if (nfull % 10000 == 0) vrb.bullet("VCF/BCF parsing [nfull=" + stb.str(nfull) + ", ncomm=" + stb.str(ncomm) + ", nrare=" + stb.str(nrare) + "]");

	}

	//Closing stuffs
	free(vgt);
	free(vac);
	free(van);
	free(vsk);
	bcf_sr_destroy(sr);
	fp_rare_bin.close();
	if (hts_close(fp_rare_bcf)) vrb.error("Non zero status when closing VCF/BCF file descriptor for rare variants");
	if (hts_close(fp_comm_bcf)) vrb.error("Non zero status when closing VCF/BCF file descriptor for common variants");


	// Report
	vrb.bullet("VCF/BCF parsing done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

