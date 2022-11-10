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


void haplotype_writer::writeHaplotypes(string fname, bool output_buffer) {
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

		if (output_buffer || (V.vec_full[vt]->bp >= input_start && V.vec_full[vt]->bp < input_stop)) {

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
					probabilities[G.GRvar_genotypes[vr][i].idx] = roundf(G.GRvar_genotypes[vr][i].prob * 1000.0) / 1000.0;
					count_alt -= 2 * major_allele;
					count_alt += a0+a1;
				}
				bcf_update_format_float(hdr, rec, "PP", probabilities, bcf_hdr_nsamples(hdr)*1);
			} else {
				for (int i = 0 ; i < H.n_samples ; i++) {
					bool a0 = H.Hvar.get(vs, 2*i+0);
					bool a1 = H.Hvar.get(vs, 2*i+1);
					count_alt += a0+a1;
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
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(G.n_samples) + " / L=" + stb.str(V.sizeFull()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(G.n_samples) + " / L=" + stb.str(V.sizeFull()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(G.n_samples) + " / L=" + stb.str(V.sizeFull()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}
