#include <containers/genotype_set.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

void genotype_set::writeValidation(std::string fname) {
	vrb.title("Write validation data in [" + fname + "]");
	tac.clock();
	std::string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	// Create VCF header
	bcf_hdr_append(hdr, std::string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr, "##source=simulVal");
	bcf_hdr_append(hdr, std::string("##contig=<ID="+ H.CHR[0] + ">").c_str());
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele ALT count\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele total count\">");


	//Add samples
	for (int i = 0 ; i < SID.size() ; i ++) bcf_hdr_add_sample(hdr, SID[i].c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	bcf_hdr_write(fp, hdr);

	//Add records
	int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
	for (int l = 0 ; l < H.n_var ; l ++) {
		bcf_clear1(rec);
		rec->rid = bcf_hdr_name2id(hdr, H.CHR[l].c_str());
		rec->pos = H.POS[l] - 1;
		bcf_update_id(hdr, rec, H.VID[l].c_str());
		std::string alleles = H.REF[l] + "," + H.ALT[l];
		bcf_update_alleles_str(hdr, rec, alleles.c_str());
		int32_t count_alt = 0, count_tot = 0;
		for (int i = 0 ; i < SID.size() ; i++) {
			genotypes[2*i+0] = bcf_gt_phased(GEN[2*i+0][l]);
			genotypes[2*i+1] = bcf_gt_phased(GEN[2*i+1][l]);
			count_alt += GEN[2*i+0][l] + GEN[2*i+1][l];
			count_tot += 2;
		}
		bcf_update_info_int32(hdr, rec, "AC", &count_alt, 1);
		bcf_update_info_int32(hdr, rec, "AN", &count_tot, 1);
		bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2);
		bcf_write1(fp, hdr, rec);
		vrb.progress("  * VCF writing", (l+1)*1.0/H.n_var);
	}
	free(genotypes);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(SID.size()) + " / L=" + stb.str(H.n_var) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(SID.size()) + " / L=" + stb.str(H.n_var) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(SID.size()) + " / L=" + stb.str(H.n_var) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}

void genotype_set::writeEstimation(std::string fname) {
	vrb.title("Write estimation data in [" + fname + "]");
	tac.clock();
	std::string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	// Create VCF header
	bcf_hdr_append(hdr, std::string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr, "##source=simulEst");
	bcf_hdr_append(hdr, std::string("##contig=<ID="+ H.CHR[0] + ">").c_str());
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele ALT count\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele total count\">");


	//Add samples
	for (int i = 0 ; i < SID.size() ; i ++) bcf_hdr_add_sample(hdr, SID[i].c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	bcf_hdr_write(fp, hdr);

	//Add records
	int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
	for (int l = 0 ; l < H.n_var ; l ++) {
		bcf_clear1(rec);
		rec->rid = bcf_hdr_name2id(hdr, H.CHR[l].c_str());
		rec->pos = H.POS[l] - 1;
		bcf_update_id(hdr, rec, H.VID[l].c_str());
		std::string alleles = H.REF[l] + "," + H.ALT[l];
		bcf_update_alleles_str(hdr, rec, alleles.c_str());
		int32_t count_alt = 0, count_tot = 0;
		for (int i = 0 ; i < SID.size() ; i++) {
			if (!MIS[i][l]) {
				if (GEN[2*i+0][l] > GEN[2*i+1][l]) {
					genotypes[2*i+0] = bcf_gt_unphased(GEN[2*i+1][l]);
					genotypes[2*i+1] = bcf_gt_unphased(GEN[2*i+0][l]);
				} else {
					genotypes[2*i+0] = bcf_gt_unphased(GEN[2*i+0][l]);
					genotypes[2*i+1] = bcf_gt_unphased(GEN[2*i+1][l]);
				}
				count_alt += GEN[2*i+0][l] + GEN[2*i+1][l];
				count_tot += 2;
			} else {
				genotypes[2*i+0] = bcf_gt_missing;
				genotypes[2*i+1] = bcf_gt_missing;
			}
		}
		bcf_update_info_int32(hdr, rec, "AC", &count_alt, 1);
		bcf_update_info_int32(hdr, rec, "AN", &count_tot, 1);
		bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2);
		bcf_write1(fp, hdr, rec);
		vrb.progress("  * VCF writing", (l+1)*1.0/H.n_var);
	}
	free(genotypes);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(SID.size()) + " / L=" + stb.str(H.n_var) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(SID.size()) + " / L=" + stb.str(H.n_var) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(SID.size()) + " / L=" + stb.str(H.n_var) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}
