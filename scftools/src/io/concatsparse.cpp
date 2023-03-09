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

#include <io/concatsparse.h>
#include <objects/rare_genotype.h>

#define BUFFERSIZE	1000

using namespace std;

concatsparse::concatsparse(string _list_txt, string _sparse_bcf, string _region, int _nthreads) {
	file_list_txt = _list_txt;
	file_rare_bcf = _sparse_bcf;
	region = _region;
	nthreads = _nthreads;
}

concatsparse::~concatsparse() {
}

void concatsparse::concat() {
	tac.clock();
	vrb.title("Concatenating multiple sparse BCF files");

	//Read files in list
	string buffer;
	vector < string > filenames;
	input_file fd_list(file_list_txt);
	if (fd_list.fail()) vrb.error("Cannot open list of VCF/BCF files");
	while (getline(fd_list, buffer)) filenames.push_back(buffer);
	fd_list.close();
	vrb.bullet(stb.str(filenames.size()) + " files detected");

	//Create filenames
	string file_rare_bin = stb.get_name_from_vcf(file_rare_bcf) + ".bin";
	string file_rare_prb = stb.get_name_from_vcf(file_rare_bcf) + ".prb";

	//Create header for output and check ordering
	int withProbs = 0;
	bcf_hdr_t * out_hdr = NULL;
	bcf1_t * line = bcf_init();
	std::vector < int > positions(filenames.size());
	std::vector < std::string > chromosomes(filenames.size());
	for (int f = 0 ; f < filenames.size() ; f ++) {
		//Open file
		htsFile *fc = hts_open(filenames[f].c_str(), "r"); if ( !fc ) vrb.error("Failed to open: " + filenames[f]);
		//Get header
		bcf_hdr_t * hdr = bcf_hdr_read(fc); if ( !hdr ) vrb.error("Failed to parse header: " + filenames[f]);
		//Merge with other headers to produce output header
		out_hdr = bcf_hdr_merge(out_hdr,hdr);
		//Check that SGEN is there
		int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "SGEN");
		if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id)) vrb.error("Cannot retrieve SGEN flag in [" + filenames[f] + "]");
		//Check is there is SPRB
		id = bcf_hdr_id2int(hdr, BCF_DT_ID, "SPRB");
		withProbs += (bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id) != 0);
		//Check chr and first position
        int ret = bcf_read(fc, hdr, line);
		if ( ret!=0 ) vrb.error("Empty file detected: " + filenames[f]);
        else { chromosomes[f] = bcf_hdr_id2name(hdr, line->rid); positions[f] = line->pos; }
		//Clean up
		bcf_hdr_destroy(hdr);
        if (hts_close(fc) != 0) vrb.error("Close failed: " + filenames[f]);
	}
    for (int i = 1; i < filenames.size() ; i++) {
    	if (chromosomes[i] != chromosomes[i-1] ) vrb.error("The files are from different chromosomes");
    	if (positions[i] < positions[i-1] ) vrb.error("The files are not in ascending order");
    }
    bool hasProbs = (withProbs == filenames.size());
    if (withProbs && !hasProbs) vrb.warning("INFO/SPRB field not found in all headers");
    if (hasProbs) vrb.bullet("Sparse probability files detected");

	//Opening sparse VCF/BCF output file
	htsFile * fp_rare_bcf = hts_open(file_rare_bcf.c_str(), "wb");
	if (!fp_rare_bcf) vrb.error("Cannot open " + file_rare_bcf + " for writing, check permissions");
	if (nthreads > 1) hts_set_threads(fp_rare_bcf, nthreads);

	//Opening Binary output file
	ofstream fp_rare_bin (file_rare_bin, std::ios::out | std::ios::binary);
	if (!fp_rare_bin) vrb.error("Cannot open " + file_rare_bin + " for writing, check permissions");

	//Opening Probs file
	ofstream fp_rare_prb;
	if (hasProbs) {
		fp_rare_prb.open(file_rare_prb, std::ios::out | std::ios::binary);
		if (!fp_rare_prb) vrb.error("Cannot open " + file_rare_prb + " for writing, check permissions");
	}

	//Write header
	bcf_hdr_add_sample(out_hdr, NULL);
	if (bcf_hdr_write(fp_rare_bcf, out_hdr)) vrb.error("Failed to write header to output BCF/VCF file");

	//Buffers for reading data
	unsigned int size_rg_buffer = BUFFERSIZE;
	unsigned int * rg_buffer = (unsigned int *)malloc(size_rg_buffer * sizeof(unsigned int));
	float * rp_buffer = (float *)malloc(size_rg_buffer * sizeof(float));
	int nsk = 0, rsk = 0, *vsk = NULL;

	//Loop over files
	ifstream fb, fp;
	uint64_t new_sgenotypes = 0;
	for (int f = 0 ; f < filenames.size() ; f ++) {
		//Verbose
		vrb.bullet("Processing [" + filenames[f] + "]");

		//Build filenames
		string frare_bcf = filenames[f];
		string frare_bin = stb.get_name_from_vcf(filenames[f]) + ".bin";
		string frare_prb = stb.get_name_from_vcf(filenames[f]) + ".prb";

		//Opening files
		htsFile * fc = hts_open(frare_bcf.c_str(), "r"); if ( !fc ) vrb.error("Failed to open: " + frare_bcf);
		fb.open(frare_bin.c_str(), std::ios::in | std::ios::binary); if (!fb) vrb.error("Failed to open: " + frare_bin);
		if (hasProbs) fp.open(frare_prb.c_str(), std::ios::in | std::ios::binary); if (!fp) vrb.error("Failed to open: " + frare_prb);

		//Read header
		bcf_hdr_t * hdr = bcf_hdr_read(fc); if ( !hdr ) vrb.error("Failed to parse header: " + frare_bcf);

		//Loop over records
		int ret;
		while (!(ret = bcf_read(fc, hdr, line))) {

			//Process sparse genotype addressing
			rsk = bcf_get_info_int32(hdr, line, "SGEN", &vsk, &nsk); if (nsk!=3) vrb.error("SGEN field is needed in rare file");
			uint64_t sgenotypes = vsk[0];
			sgenotypes *= MOD30BITS;
			sgenotypes += vsk[1];
			uint32_t ngenotypes = vsk[2];
			if (ngenotypes > size_rg_buffer) {
				size_rg_buffer = ngenotypes;
				rg_buffer = (unsigned int *)realloc(rg_buffer, size_rg_buffer * sizeof(unsigned int));
				rp_buffer = (float *)realloc(rp_buffer, size_rg_buffer * sizeof(float));
			}

			//Read sparse genotypes in file
			fb.seekg(sgenotypes * sizeof(unsigned int), fb.beg);
			fb.read(reinterpret_cast < char * > (rg_buffer), ngenotypes * sizeof(unsigned int));

			//Read sparse probabilities
			if (hasProbs) {
				fp.seekg(sgenotypes * sizeof(unsigned int), fp.beg);
				fp.read(reinterpret_cast < char * > (rp_buffer), ngenotypes * sizeof(float));
			}

			//Copy over sparse data from src to dest
			fp_rare_bin.write(reinterpret_cast < char * > (rg_buffer), ngenotypes * sizeof(unsigned int));

			//Copy over sparse probs from src to dest
			if (hasProbs) fp_rare_prb.write(reinterpret_cast < char * > (rp_buffer), ngenotypes * sizeof(float));

			//Update BCF record with new data
			vsk[0] = new_sgenotypes / MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
			vsk[1] = new_sgenotypes % MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
			bcf_update_info_int32(out_hdr, line, "SGEN", vsk, 3);
			if (hasProbs) bcf_update_info_int32(out_hdr, line, "SPRB", vsk, 3);
			new_sgenotypes += ngenotypes;

			//Write BCF record in output file
			if (bcf_write1(fp_rare_bcf, out_hdr, line) < 0) vrb.error("Failing to write VCF/record");
		}

		//Close
		bcf_hdr_destroy(hdr);
		if (hts_close(fc) != 0) vrb.error("Close failed: " + filenames[f]);
		fb.close();
		fp.close();

	}


	//Closing stuffs
	free(rg_buffer);
	free(rp_buffer);
	fp_rare_bin.close();
	fp_rare_prb.close();
	if (hts_close(fp_rare_bcf)) vrb.error("Non zero status when closing VCF/BCF file descriptor for rare variants");

	// Report
	vrb.bullet("Sparse VCF/BCF concat done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	//Indexing
	vrb.bullet("Indexing [" + file_rare_bcf + "]");
	if (bcf_index_build3(file_rare_bcf.c_str(), NULL, 14, nthreads) < 0) vrb.error("Fail to index file");
}
