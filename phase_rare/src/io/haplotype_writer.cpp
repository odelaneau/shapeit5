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
#include <io/haplotype_writer.h>

using namespace std;

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

haplotype_writer::haplotype_writer(haplotype_set & _H, genotype_set & _G, variant_map & _V, int _nthreads): H(_H), G(_G), V(_V) {
	nthreads = _nthreads;
}

haplotype_writer::~haplotype_writer() {
}


void haplotype_writer::setRegions(int _input_start, int _input_stop) {
	input_start = _input_start;
	input_stop = _input_stop;
}

void haplotype_writer::writeHaplotypesPlain(string fname, bool output_buffer) {
	// Init
	tac.clock();
	string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	if (nthreads > 1) hts_set_threads(fp, nthreads);
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	// Create VCF header
	bcf_hdr_append(hdr, string("##source=shapeit5 phase 2 v" + string(PHASE2_VERSION)).c_str());
	bcf_hdr_append(hdr, string("##contig=<ID="+ V.vec_full[0]->chr + ">").c_str());
	bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=PP,Number=1,Type=Float,Description=\"Phasing confidence\">");

	//Add samples
	for (int i = 0 ; i < G.n_samples ; i ++) bcf_hdr_add_sample(hdr, G.names[i].c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	if (bcf_hdr_write(fp, hdr) < 0) vrb.error("Failing to write VCF/header");

	//Add records
	int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
	float * probabilities = (float*)malloc(bcf_hdr_nsamples(hdr)*1*sizeof(float));

	for (int vt = 0, vc = 0, vs = 0, vr = 0 ; vt < V.sizeFull() ; vt ++) {

		if (output_buffer || (V.vec_full[vt]->bp >= input_start && V.vec_full[vt]->bp <= input_stop)) {

			//Variant informations
			bcf_clear1(rec);
			rec->rid = bcf_hdr_name2id(hdr, V.vec_full[vt]->chr.c_str());
			rec->pos = V.vec_full[vt]->bp - 1;
			bcf_update_id(hdr, rec, V.vec_full[vt]->id.c_str());
			string alleles = V.vec_full[vt]->ref + "," + V.vec_full[vt]->alt;
			bcf_update_alleles_str(hdr, rec, alleles.c_str());

			//Genotypes
			int count_alt = 0;

			if (V.vec_full[vt]->type == VARTYPE_RARE) {
				bool major_allele = !V.vec_full[vt]->minor;
				for (int i = 0 ; i < G.n_samples ; i++) {
					genotypes[2*i+0] = bcf_gt_phased(major_allele);
					genotypes[2*i+1] = bcf_gt_phased(major_allele);
					bcf_float_set_missing(probabilities[i]);
					count_alt += 2 * major_allele;
				}
				for (int i = 0 ; i < G.GRvar_genotypes[vr].size() ; i++) {
					bool a0 = G.GRvar_genotypes[vr][i].al0;
					bool a1 = G.GRvar_genotypes[vr][i].al1;
					genotypes[2*G.GRvar_genotypes[vr][i].idx+0] = bcf_gt_phased(a0);
					genotypes[2*G.GRvar_genotypes[vr][i].idx+1] = bcf_gt_phased(a1);
					probabilities[G.GRvar_genotypes[vr][i].idx] = roundf(min(G.GRvar_genotypes[vr][i].prob, 1.0f) * 1000.0) / 1000.0;
					count_alt -= 2 * major_allele;
					count_alt += a0+a1;
				}
				bcf_update_format_float(hdr, rec, "PP", probabilities, bcf_hdr_nsamples(hdr)*1);
			} else {
				for (int i = 0 ; i < H.n_samples ; i++) {
					bool a0 = H.Hvar.get(vs, 2*i+0);
					bool a1 = H.Hvar.get(vs, 2*i+1);
					genotypes[2*i+0] = bcf_gt_phased(a0);
					genotypes[2*i+1] = bcf_gt_phased(a1);
					count_alt += a0+a1;
				}
			}

			bcf_update_info_int32(hdr, rec, "AC", &count_alt, 1);
			bcf_update_info_int32(hdr, rec, "AN", &G.n_samples, 1);
			bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2);
			if (V.vec_full[vt]->type == VARTYPE_RARE) bcf_update_format_float(hdr, rec, "PP", probabilities, bcf_hdr_nsamples(hdr)*1);

			if (bcf_write1(fp, hdr, rec) < 0) vrb.error("Failing to write VCF/record");
		}

		switch (V.vec_full[vt]->type) {
		case VARTYPE_SCAF :	vs++; break;
		case VARTYPE_COMM :	vc++; break;
		case VARTYPE_RARE :	vr++; break;
		}

		vrb.progress("  * VCF writing", (vt+1)*1.0/V.sizeFull());
	}
	free(genotypes);
	free(probabilities);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(G.n_samples) + " / L=" + stb.str(V.sizeFull()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(G.n_samples) + " / L=" + stb.str(V.sizeFull()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(G.n_samples) + " / L=" + stb.str(V.sizeFull()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}

	vrb.bullet("Indexing ["+fname + "]");
	if (bcf_index_build3(fname.c_str(), NULL, 14, nthreads) < 0) vrb.error("Fail to index file");
}

void haplotype_writer::writeHaplotypesSparse(string fname) {

	// Init
	tac.clock();

	//Opening all sparse files
	string file_rare_bcf = fname;
	string file_rare_bin = stb.get_name_from_vcf(fname) + ".bin";
	string file_rare_prb = stb.get_name_from_vcf(fname) + ".prb";
	htsFile * fp_rare_bcf = hts_open(file_rare_bcf.c_str(), "wb");
	if (!fp_rare_bcf) vrb.error("Cannot open " + file_rare_bcf + " for writing, check permissions");
	ofstream fp_rare_bin (file_rare_bin, std::ios::out | std::ios::binary);
	if (!fp_rare_bin) vrb.error("Cannot open " + file_rare_bin + " for writing, check permissions");
	ofstream fp_rare_prb (file_rare_prb, std::ios::out | std::ios::binary);
	if (!fp_rare_prb) vrb.error("Cannot open " + file_rare_bin + " for writing, check permissions");

	//Create and write VCF header for rare BCF
	bcf_hdr_t * hdr_rare_bcf = bcf_hdr_init("w");
	bcf_hdr_append(hdr_rare_bcf, string("##source=shapeit5 convert v" + string(SCFTLS_VERSION)).c_str());
	bcf_hdr_append(hdr_rare_bcf, string("##contig=<ID="+ V.vec_full[0]->chr + ">").c_str());
	bcf_hdr_append(hdr_rare_bcf, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr_rare_bcf, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr_rare_bcf, "##INFO=<ID=SGEN,Number=3,Type=Integer,Description=\"Index and number of sparse genotypes in sparse file\">");
	bcf_hdr_append(hdr_rare_bcf, "##INFO=<ID=SPRB,Number=3,Type=Integer,Description=\"Index and number of sparse probabilities in sparse file\">");
	bcf_hdr_add_sample(hdr_rare_bcf, NULL);
	if (bcf_hdr_write(fp_rare_bcf, hdr_rare_bcf) < 0) vrb.error("Failing to write VCF/header for rare variants");

	//
	bcf1_t *rec = bcf_init1();
	uint64_t seek = 0;
	int nsk, rsk, *vsk = (int*)malloc(3*sizeof(int));
	unsigned int * genotypes = (unsigned int *)malloc(G.n_samples*sizeof(unsigned int));
	float * probabilities = (float*)malloc(G.n_samples*sizeof(float));

	for (int vr = 0 ; vr < V.sizeRare() ; vr ++) {

		//Variant informations
		bcf_clear1(rec);
		rec->rid = bcf_hdr_name2id(hdr_rare_bcf, V.vec_rare[vr]->chr.c_str());
		rec->pos = V.vec_rare[vr]->bp - 1;
		bcf_update_id(hdr_rare_bcf, rec, V.vec_rare[vr]->id.c_str());
		string alleles = V.vec_rare[vr]->ref + "," + V.vec_rare[vr]->alt;
		bcf_update_alleles_str(hdr_rare_bcf, rec, alleles.c_str());

		//Seek information
		vsk[0] = seek / MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
		vsk[1] = seek % MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
		vsk[2] = 0;

		//Compact genotypes and probabilities
		bool major_allele = !V.vec_rare[vr]->minor;
		unsigned int count_alt = 2 * G.n_samples * major_allele;
		for (int r = 0 ; r < G.GRvar_genotypes[vr].size() ; r ++) {
			bool a0 = G.GRvar_genotypes[vr][r].al0;
			bool a1 = G.GRvar_genotypes[vr][r].al1;
			if (a0 != major_allele || a1 != major_allele) {
				probabilities[r] = min(G.GRvar_genotypes[vr][r].prob, 1.0f);
				rare_genotype rg_struct = rare_genotype(G.GRvar_genotypes[vr][r].idx, (a0!=a1), 0, a0, a1, 1);
				genotypes[r] = rg_struct.get();
				count_alt -= 2 * major_allele;
				count_alt += a0+a1;
				vsk[2] ++;
				seek ++;
			}
		}

		//Write genotypes and probabilities
		fp_rare_bin.write(reinterpret_cast < char * > (genotypes), vsk[2] * sizeof(unsigned int));
		fp_rare_prb.write(reinterpret_cast < char * > (probabilities), vsk[2] * sizeof(float));

		//BCF INFO fields
		bcf_update_info_int32(hdr_rare_bcf, rec, "AC", &count_alt, 1);
		bcf_update_info_int32(hdr_rare_bcf, rec, "AN", &G.n_samples, 1);
		bcf_update_info_int32(hdr_rare_bcf, rec, "SGEN", vsk, 3);
		bcf_update_info_int32(hdr_rare_bcf, rec, "SPRB", vsk, 3);
		if (bcf_write1(fp_rare_bcf, hdr_rare_bcf, rec) < 0) vrb.error("Failing to write VCF/record");

		vrb.progress("  * VCF writing", (vr+1)*1.0/V.sizeRare());
	}
	free(vsk);
	free(probabilities);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr_rare_bcf);
	if (hts_close(fp_rare_bcf)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	fp_rare_bin.close();
	fp_rare_prb.close();
	vrb.bullet("Sparse BCF writing [N=" + stb.str(G.n_samples) + " / L=" + stb.str(V.sizeRare()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)");

	vrb.bullet("Indexing ["+file_rare_bcf + "]");
	if (bcf_index_build3(file_rare_bcf.c_str(), NULL, 14, nthreads) < 0) vrb.error("Fail to index file");
}

