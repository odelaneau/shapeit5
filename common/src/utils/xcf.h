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

#ifndef _XCF_H
#define _XCF_H

#include <utils/otools.h>

#define FILE_VOID	0	//No data
#define FILE_BCF	1	//Data in BCF file
#define FILE_BINARY	2	//Data in Binary file

#define RECORD_VOID				0
#define RECORD_BCFVCF_GENOTYPE	1
#define RECORD_SPARSE_GENOTYPE	2
#define RECORD_SPARSE_HAPLOTYPE	3
#define RECORD_BINARY_GENOTYPE	4
#define RECORD_BINARY_HAPLOTYPE	5
#define RECORD_NUMBER_TYPES		6
/*****************************************************************************/
/*****************************************************************************/
/******						XCF_READER									******/
/*****************************************************************************/
/*****************************************************************************/

class xcf_reader {
public:
	//HTS part
	uint32_t sync_number;
	bcf_srs_t * sync_reader;
	std::vector < bcf1_t * > sync_lines;
	std::vector < int32_t > sync_types;			//Type of data: [DATA_EMPTY, FILE_BCF, FILE_BINARY]
	std::vector < bool > sync_flags;			//Has record?

	//Variant information
	std::string chr;
	uint32_t pos;
	std::string ref;
	std::string alt;
	std::string rsid;
	std::vector < uint32_t > AC;
	std::vector < uint32_t > AN;

	//INFO field
	int32_t nAC, nAN, nSK;
	int32_t * vAC;
	int32_t * vAN;
	int32_t * vSK;

	//Pedigree file
	std::vector < uint32_t> ind_number;
	std::vector < std::vector < std::string > > ind_names;
	std::vector < std::vector < std::string > > ind_fathers;
	std::vector < std::vector < std::string > > ind_mothers;

	//Binary files [files x types]
	std::vector < std::ifstream > bin_fds;		//File Descriptors
	std::vector < int32_t > bin_type;			//Type of Binary record					//Integer 1 in INFO/SEEK field
	std::vector < uint64_t > bin_seek;			//Location of Binary record				//Integer 2 and 3 in INFO/SEEK field
	std::vector < uint32_t > bin_size;			//Amount of Binary records in bytes		//Integer 4 in INFO/SEEK field

	//CONSTRUCTOR
	xcf_reader(std::string region, uint32_t nthreads) {
		sync_number = 0;
		sync_reader = bcf_sr_init();
		sync_reader->collapse = COLLAPSE_NONE;
		sync_reader->require_index = 1;
		if (nthreads > 1) bcf_sr_set_threads(sync_reader, nthreads);
		if (bcf_sr_set_regions(sync_reader, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "]");
		vAC = vAN = vSK = NULL;
		nAC = nAN = nSK = 0;
	}

	//DESTRUCTOR
	~xcf_reader() {
		//bcf_sr_destroy(sync_reader);
	}

	//ADD A NEW FILE IN THE SYNCHRONIZED READER
	int32_t addFile(std::string fname) {
		std::string buffer;
		std::vector < std::string > tokens;

		//Open BCF file and add it in the synchronized reader
		if (!(bcf_sr_add_reader (sync_reader, fname.c_str()))) {
			switch (sync_reader->errnum) {
			case not_bgzf:			vrb.error("Opening [" + fname + "]: not compressed with bgzip"); break;
			case idx_load_failed: 	vrb.error("Opening [" + fname + "]: impossible to load index file"); break;
			case file_type_error: 	vrb.error("Opening [" + fname + "]: file format not supported by HTSlib"); break;
			default : 				vrb.error("Opening [" + fname + "]: unknown error"); break;
			}
		}

		//Allocation for Binary record information
		sync_lines.push_back(bcf_init());
		sync_types.push_back(FILE_VOID);
		sync_flags.push_back(false);
		bin_fds.push_back(std::ifstream());
		bin_type.push_back(0);
		bin_seek.push_back(0);
		bin_size.push_back(0);
		AC.push_back(0);
		AN.push_back(0);

		//Check header for associated binary file
		int32_t flagSEEK = bcf_hdr_idinfo_exists(sync_reader->readers[sync_number].header, BCF_HL_INFO, bcf_hdr_id2int(sync_reader->readers[sync_number].header, BCF_DT_ID, "SEEK"));

		//Check header for number of samples in BCF
		uint32_t nsamples = bcf_hdr_nsamples(sync_reader->readers[sync_number].header);

		/************************************************************************************/
		/*   CASE1: There is a binary file and no data in the BCF / Open the binary file	*/
		/************************************************************************************/
		if (flagSEEK && nsamples == 0) {
			//Open Binary file
			std::string bfname = stb.get_name_from_vcf(fname) + ".bin";
			bin_fds[sync_number].open(bfname.c_str(), std::ios::in | std::ios::binary);
			if (!bin_fds[sync_number]) vrb.error("Cannot open file [" + bfname + "] for reading");
			//Read PED file
			std::string ped_fname = stb.get_name_from_vcf(fname) + ".fam";
			input_file fdp(ped_fname);
			if (fdp.fail()) vrb.error("Cannot open pedigree file [" + ped_fname + "] for reading");
			ind_names.push_back(std::vector < std::string >());
			ind_fathers.push_back(std::vector < std::string >());
			ind_mothers.push_back(std::vector < std::string >());
			while (getline(fdp, buffer)) {
				stb.split(buffer, tokens);
				ind_names[sync_number].push_back(tokens[0]);
				if (tokens.size() >=3 ) { ind_fathers[sync_number].push_back(tokens[1]); ind_mothers[sync_number].push_back(tokens[2]); }
				else { ind_fathers[sync_number].push_back("NA"); ind_mothers[sync_number].push_back("NA"); }
			}
			ind_number.push_back(ind_names[sync_number].size());
			sync_types[sync_number] = FILE_BINARY;
			fdp.close();
		}

		/************************************************************************************/
		/*   CASE2: There is NOT a binary file and data in the BCF 							*/
		/************************************************************************************/
		if (!flagSEEK && nsamples != 0) {
			ind_names.push_back(std::vector < std::string >());
			for (int i = 0 ; i < nsamples ; i ++)
				ind_names[sync_number].push_back(std::string(sync_reader->readers[sync_number].header->samples[i]));
			ind_number.push_back(ind_names[sync_number].size());
			ind_fathers.push_back(std::vector < std::string >(ind_number[sync_number], "NA"));
			ind_mothers.push_back(std::vector < std::string >(ind_number[sync_number], "NA"));
			sync_types[sync_number] = FILE_BCF;
		}

		/************************************************************************************/
		/*   CASE3: There is NOT a binary file and No data in the BCF 						*/
		/************************************************************************************/
		if (!flagSEEK && nsamples == 0) {
			ind_number.push_back(0);
			sync_types[sync_number] = FILE_VOID;
		}

		/************************************************************************************/
		/*   CASE4: There is a binary file and data in the BCF 								*/
		/************************************************************************************/
		if (flagSEEK && nsamples != 0) vrb.error("Binary file found for a non-empty BCF file");

		//Increment number of readers
		sync_number++;
		return (sync_number-1);
	}

	//GET SAMPLES
	int32_t getSamples(uint32_t file, std::vector < std::string > & samples) {
		samples.clear();
		for (uint32_t i = 0 ; i < ind_number[file] ; i ++) samples.push_back(ind_names[file][i]);
		return samples.size();
	}

	//GET ACs/ANs/AFs
	uint32_t getAC(uint32_t file) { return AC[file]; }
	uint32_t getAN(uint32_t file) { return AN[file]; }
	uint32_t getAC() { return std::accumulate(AC.begin(), AC.end(), 0); }
	uint32_t getAN() { return std::accumulate(AN.begin(), AN.end(), 0); }
	float getAF(uint32_t file) { return AC[file]*1.0f/AN[file]; }
	float getAF() { return std::accumulate(AC.begin(), AC.end(), 0)*1.0f/std::accumulate(AN.begin(), AN.end(), 0); }

	//SET SYNCHRONIZED READER TO NEXT RECORD
	int32_t nextRecord() {

	next_record:
		//Go to next record
		int32_t ret = bcf_sr_next_line (sync_reader);
		if (!ret) return 0;

		//Initialize
		std::fill(sync_flags.begin(), sync_flags.end(), false);
		std::fill(AC.begin(), AC.end(), 0);
		std::fill(AN.begin(), AN.end(), 0);
		std::fill(bin_type.begin(), bin_type.end(), RECORD_VOID);
		std::fill(bin_seek.begin(), bin_seek.end(), 0);
		std::fill(bin_size.begin(), bin_size.end(), 0);

		//Loop over readers
		for (uint32_t r = 0, firstfile = 0 ; r < sync_number ; r++) {

			//Check if reader has a record
			bool hasRecord = bcf_sr_has_line(sync_reader, r);

			//First time we see a record for this variant
			if (hasRecord) {

				//Get the record
				sync_lines[sync_number] = bcf_sr_get_line(sync_reader, 0);

				//If non bi-allelic, goto next variant
				if (sync_lines[sync_number]->n_allele != 2) goto next_record;

				//Get variant information
				if (!firstfile) {
					chr = bcf_hdr_id2name(sync_reader->readers[0].header, sync_lines[sync_number]->rid);
					pos = sync_lines[sync_number]->pos + 1;
					rsid = std::string(sync_lines[sync_number]->d.id);
					ref = std::string(sync_lines[sync_number]->d.allele[0]);
					alt = std::string(sync_lines[sync_number]->d.allele[1]);
					firstfile = 1;
				}

				//Get AC/AN information
				int32_t rAC = bcf_get_info_int32(sync_reader->readers[r].header, sync_lines[sync_number], "AC", &vAC, &nAC);
				int32_t rAN = bcf_get_info_int32(sync_reader->readers[r].header, sync_lines[sync_number], "AN", &vAN, &nAN);
				if (rAC != 1) vrb.error("AC field is needed in main file for MAF filtering");
				if (rAN != 1) vrb.error("AN field is needed in main file for MAF filtering");
				AC[r] = vAC[0]; AN[r] = vAN[0];

				//Get SEEK information
				if (sync_types[r] == FILE_BINARY) {
					int32_t rsk = bcf_get_info_int32(sync_reader->readers[r].header, sync_lines[sync_number], "SEEK", &vSK, &nSK);
					if (nSK != 4) vrb.error("INFO/SEEK field should contain 4 numbers");
					else {
						bin_type[r] = vSK[0];
						bin_seek[r] = vSK[1];
						bin_seek[r] *= MOD30BITS;
						bin_seek[r] += vSK[2];
						bin_size[r] = vSK[3];
					}
				} else if (sync_types[r] == FILE_BCF) {
					bin_type[r] = RECORD_BCFVCF_GENOTYPE;
					bin_seek[r] = 0;
					bin_size[r] = 0;
				}

				//Set flag
				sync_flags[r] = true;
			}
		}

		//Return number of files with a record
		return ret;
	}

	//CHECK IF FILE HAS A RECORD THERE
	int32_t hasRecord(uint32_t file) {
		return (sync_flags[file]);
	}

	//CHECK THE TYPE OF FILE
	// 0: BCF file without sample data available
	// 1: BCF file with sample data available
	// 2 []: Binary file There is a record with sample data in Binary file / return the type of binary record
	int32_t typeFile(uint32_t file) {
		return sync_types[file];
	}

	//CHECK THE TYPE OF RECORD
	// 0: BCF file without sample data available
	// 1: BCF file with sample data available
	// 2 []: Binary file There is a record with sample data in Binary file / return the type of binary record
	int32_t typeRecord(uint32_t file) {
		return bin_type[file];
	}

	//READ DATA OF THE AVAILABLE RECORD
	// =0: No sample data available
	// >0: Amount of data read in bytes
	int32_t readRecord(uint32_t file, char ** buffer) {

		//No data in this file for the current record
		if (!sync_flags[file]) return 0;

		//No data since file is empty
		if (sync_types[file] == FILE_VOID) return 0;

		//Data is in BCF file
		if (sync_types[file] == FILE_BCF) {
			//Read genotypes [assuming buffer to be allocated!]
			int32_t ndp = ind_number[file]*2;
			int32_t rdp = bcf_get_genotypes(sync_reader->readers[file].header, sync_lines[sync_number], buffer, &ndp);
			assert(rdp == (ind_number[file]*2));
			return ndp * sizeof(int32_t);
		}

		//Data is in binary file
		else {
			//Seek to right position
			bin_fds[file].seekg(bin_seek[file], bin_fds[file].beg);

			//Read data in Binary file
			bin_fds[file].read(*buffer, bin_size[file]);

			//Return amount of data in bytes read in file
			return bin_size[file];
		}
	}

	void close() {
		free(vSK); free(vAC); free(vAN);
		for (uint32_t r = 0 ; r < sync_number ; r++) if (sync_types[r]>=2) bin_fds[r].close();
		bcf_sr_destroy(sync_reader);
	}
};


/*****************************************************************************/
/*****************************************************************************/
/******						XCF_WRITER									******/
/*****************************************************************************/
/*****************************************************************************/

class xcf_writer {
public:
	//HTS part
	std::string hts_fname;
	htsFile * hts_fd;
	bcf_hdr_t * hts_hdr;
	bcf1_t * hts_record;
	bool hts_genotypes;
	uint32_t nthreads;

	//INFO field
	int32_t nsk;
	int32_t rsk;
	int32_t * vsk;

	//Pedigree file
	uint32_t ind_number;
	std::vector < std::string > ind_names;
	std::vector < std::string > ind_fathers;
	std::vector < std::string > ind_mothers;

	//Binary files [files x types]
	std::ofstream bin_fds;						//File Descriptors
	uint32_t bin_type;							//Type of Binary record					//Integer 1 in INFO/SEEK field
	uint64_t bin_seek;							//Location of Binary record				//Integer 2 and 3 in INFO/SEEK field
	uint32_t bin_size;							//Amount of Binary records in bytes		//Integer 4 in INFO/SEEK field

	//CONSTRUCTOR
	xcf_writer(std::string fname, bool _hts_genotypes, uint32_t _nthreads) {
		std::string oformat;
		hts_fname = fname;

		if (hts_fname == "-") {
			oformat = "wbu";		//Uncompressed BCF for stdout
		} else {
			oformat = "wb";			//Compressed BCF for file
			if (hts_fname.size() < 5) vrb.error("Filename too short");
			if (hts_fname.substr(fname.size()-3) != "bcf") vrb.error("Filename extension should be *.bcf");
		}

		hts_genotypes = _hts_genotypes;
		nthreads = _nthreads;
		bin_type = 0;
		bin_seek = 0;
		bin_size = 0;
		hts_record = bcf_init1();
		vsk = (int32_t *)malloc(4 * sizeof(int32_t *));
		nsk = rsk = 0;

		hts_fd = hts_open(hts_fname.c_str(), oformat.c_str());
		if (nthreads > 1) hts_set_threads(hts_fd, nthreads);

		if (!hts_genotypes) {
			//BINARY
			std::string bfname = stb.get_name_from_vcf(hts_fname) + ".bin";
			bin_fds.open(bfname.c_str(), std::ios::out | std::ios::binary);
			if (!bin_fds) vrb.error("Cannot open file [" + bfname + "] for writing");
		}
	}

	//DESTRUCTOR
	~xcf_writer() {
		//close();
	}

	//Write sample IDs
	void writeHeader(std::vector < std::string > & samples, std::string contig) {
		//
		hts_hdr = bcf_hdr_init("w");
		bcf_hdr_append(hts_hdr, std::string("##fileDate="+tac.date()).c_str());
		bcf_hdr_append(hts_hdr, std::string("##source=XCFtools").c_str());
		bcf_hdr_append(hts_hdr, std::string("##contig=<ID="+ contig + ">").c_str());
		bcf_hdr_append(hts_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ALT allele count\">");
		bcf_hdr_append(hts_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of alleles\">");
		if (hts_genotypes) bcf_hdr_append(hts_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
		else bcf_hdr_append(hts_hdr, "##INFO=<ID=SEEK,Number=4,Type=Integer,Description=\"SEEK binary file information\">");

		//Write sample IDs
		if (hts_genotypes) {
			//Samples are in BCF header
			for (uint32_t i = 0 ; i < samples.size() ; i++) bcf_hdr_add_sample(hts_hdr, samples[i].c_str());
			bcf_hdr_add_sample(hts_hdr, NULL);      // to update internal structures
		} else {
			//Samples are in PED file
			std::string ffname = stb.get_name_from_vcf(hts_fname) + ".fam";
			output_file fd (ffname);
			if (fd.fail()) vrb.error("Cannot open [" + ffname + "] for writing");
			for (uint32_t i = 0 ; i < samples.size() ; i++) fd << samples[i] << "\tNA\tNA" << std::endl;
			fd.close();
		}
		if (bcf_hdr_write(hts_fd, hts_hdr) < 0) vrb.error("Failing to write VCF/header");
	}

	//Write sample IDs
	void writeHeader(std::vector < std::string > & samples, std::vector < int > & fathers, std::vector < int > & mothers, std::string contig) {
		//
		hts_hdr = bcf_hdr_init("w");
		bcf_hdr_append(hts_hdr, std::string("##fileDate="+tac.date()).c_str());
		bcf_hdr_append(hts_hdr, std::string("##source=XCFtools").c_str());
		bcf_hdr_append(hts_hdr, std::string("##contig=<ID="+ contig + ">").c_str());
		bcf_hdr_append(hts_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ALT allele count\">");
		bcf_hdr_append(hts_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of alleles\">");
		if (hts_genotypes) bcf_hdr_append(hts_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
		else bcf_hdr_append(hts_hdr, "##INFO=<ID=SEEK,Number=4,Type=Integer,Description=\"SEEK binary file information\">");

		//Write sample IDs
		if (hts_genotypes) {
			//Samples are in BCF header
			for (uint32_t i = 0 ; i < samples.size() ; i++) bcf_hdr_add_sample(hts_hdr, samples[i].c_str());
			bcf_hdr_add_sample(hts_hdr, NULL);      // to update internal structures
		} else {
			//Samples are in PED file
			std::string ffname = stb.get_name_from_vcf(hts_fname) + ".fam";
			output_file fd (ffname);
			if (fd.fail()) vrb.error("Cannot open [" + ffname + "] for writing");
			for (uint32_t i = 0 ; i < samples.size() ; i++) {
				std::string sfather = (fathers[i]>=0)?samples[fathers[i]]:"NA";
				std::string smother = (mothers[i]>=0)?samples[mothers[i]]:"NA";
				fd << samples[i] << "\t" << sfather << "\t" << smother << std::endl;
			}
			fd.close();
		}
		if (bcf_hdr_write(hts_fd, hts_hdr) < 0) vrb.error("Failing to write VCF/header");
	}

	//Write variant information
	void writeInfo(std::string chr, uint32_t pos, std::string ref, std::string alt, std::string rsid, uint32_t AC, uint32_t AN) {
		bcf_clear1(hts_record);
		hts_record->rid = bcf_hdr_name2id(hts_hdr, chr.c_str());
		hts_record->pos = pos - 1;
		bcf_update_id(hts_hdr, hts_record, rsid.c_str());
		std::string alleles = ref + "," + alt;
		bcf_update_alleles_str(hts_hdr, hts_record, alleles.c_str());
		bcf_update_info_int32(hts_hdr, hts_record, "AC", &AC, 1);
		bcf_update_info_int32(hts_hdr, hts_record, "AN", &AN, 1);
	}

	//Write genotypes
	void writeRecord(uint32_t type, char * buffer, uint32_t nbytes) {
		if (hts_genotypes) {
			bcf_update_genotypes(hts_hdr, hts_record, buffer, nbytes/sizeof(int32_t));
		} else {
			vsk[0] = type;
			vsk[1] = bin_seek / MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
			vsk[2] = bin_seek % MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
			vsk[3] = nbytes;
			bin_fds.write(buffer, nbytes);
			bin_seek += nbytes;
			bcf_update_info_int32(hts_hdr, hts_record, "SEEK", vsk, 4);
		}
		if (bcf_write1(hts_fd, hts_hdr, hts_record) < 0) vrb.error("Failing to write VCF/record for rare variants");
	}

	void close() {
		free(vsk);
		bcf_destroy1(hts_record);
		bcf_hdr_destroy(hts_hdr);
		if (hts_close(hts_fd)) vrb.error("Non zero status when closing [" + hts_fname + "]");
		if (hts_fname!="-") {
			if (bcf_index_build3(hts_fname.c_str(), NULL, 14, nthreads) < 0) vrb.error("Fail to index file [" + hts_fname + "]");
		}

	}
};

#endif

