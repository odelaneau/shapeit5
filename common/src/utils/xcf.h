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

#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>
#include <regex>
#include <chrono>
#include <iomanip>


//INCLUDE HTS LIBRARY
extern "C" {
	#include <htslib/synced_bcf_reader.h>
}

#define FILE_VOID	0					//No data
#define FILE_BCF	1					//Data in BCF file
#define FILE_BINARY	2					//Data in Binary file

#define RECORD_VOID				0		//No record
#define RECORD_BCFVCF_GENOTYPE	1		//Record in BCF GT format
#define RECORD_SPARSE_GENOTYPE	2		//Record in sparse genotype format (see rare_genotype.h)
#define RECORD_SPARSE_HAPLOTYPE	3		//Record in sparse haplotype format (uint32_t for indexing)
#define RECORD_BINARY_GENOTYPE	4		//Record in binary genotype format (2bits per genotype; 10 for missing)
#define RECORD_BINARY_HAPLOTYPE	5		//Record in binary haplotype format (1bit per allele; no missing allowed)
#define RECORD_NUMBER_TYPES		6

#define MOD30BITS			0x40000000

/*****************************************************************************/
/*****************************************************************************/
/******						XCF_UTILS									******/
/*****************************************************************************/
/*****************************************************************************/

namespace helper_tools
{
	inline std::string findExtension ( const std::string & filename ) {
	   auto position = filename.find_last_of ( '.' ) ;
	   if ( position == std::string::npos )
	      return "" ;
	   else {
	      std::string extension ( filename.substr( position + 1 ) ) ;
	      if (std::regex_search (extension, std::regex("[^A-Za-z0-9]") ))
	         return "" ;
	      else
	         return extension ;
	   }
	}

	inline std::string get_name_from_vcf(std::string filename)
	{
		std::string ext = findExtension(filename);
		if (ext == "vcf" || "bcf")
		{
			size_t lastdot = filename.find_last_of(".");
			if (lastdot == std::string::npos) return filename;
			return filename.substr(0, lastdot);
		}
		else if (ext=="gz") //check for vcf.gz
		{
			size_t lastdot = filename.find_last_of(".");
			if (lastdot == std::string::npos) return filename;
			std::string filename2 =  filename.substr(0, lastdot);
			if (findExtension(filename2) == "vcf")
			{
				lastdot = filename2.find_last_of(".");
				if (lastdot == std::string::npos) return filename2;
				return filename.substr(0, lastdot);
			}
		}
		return filename;
	}

	inline int split(const std::string & str, std::vector < std::string > & tokens, std::string sep = " 	", unsigned int n_max_tokens = 1000000) {
		tokens.clear();
		if (str == ""){
			tokens.push_back("");
			return tokens.size();
		}
		std::string::size_type p_last = str.find_first_not_of(sep, 0);
		std::string::size_type p_curr = str.find_first_of(sep, p_last);
		while ((std::string::npos != p_curr || std::string::npos != p_last) && tokens.size() < n_max_tokens) {
			tokens.push_back(str.substr(p_last, p_curr - p_last));
			p_last = str.find_first_not_of(sep, p_curr);
			p_curr = str.find_first_of(sep, p_last);
		}
		if (tokens.back()[tokens.back().size()-1] == '\r') tokens.back() = tokens.back().substr(0, tokens.back().size()-1);
		return tokens.size();
	}

	inline void error(std::string s) {
		std::cout << std::endl << "\x1B[31m" << "ERROR: " <<  "\033[0m" << s << std::endl;
		exit(EXIT_FAILURE);
	}

	inline std::string date() {
		auto now = std::chrono::system_clock::now();
		auto in_time_t = std::chrono::system_clock::to_time_t(now);
		std::stringstream ss;
	    ss << std::put_time(std::localtime(&in_time_t), "%d/%m/%Y - %X");
	    return ss.str();
	}

}

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
	bool multi;
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
	xcf_reader(std::string region, uint32_t nthreads) : pos(0), multi(false) {
		sync_number = 0;
		sync_reader = bcf_sr_init();
		sync_reader->collapse = COLLAPSE_NONE;
		sync_reader->require_index = 1;
		if (nthreads > 1) bcf_sr_set_threads(sync_reader, nthreads);
		if (bcf_sr_set_regions(sync_reader, region.c_str(), 0) == -1) helper_tools::error("Impossible to jump to region [" + region + "]");
		if (bcf_sr_set_targets(sync_reader, region.c_str(), 0, 0) == -1) helper_tools::error("Impossible to constrain to region [" + region + "]");
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
			case not_bgzf:			helper_tools::error("Opening [" + fname + "]: not compressed with bgzip"); break;
			case idx_load_failed: 	helper_tools::error("Opening [" + fname + "]: impossible to load index file"); break;
			case file_type_error: 	helper_tools::error("Opening [" + fname + "]: file format not supported by HTSlib"); break;
			default : 				helper_tools::error("Opening [" + fname + "]: unknown error"); break;
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
			std::string bfname = helper_tools::get_name_from_vcf(fname) + ".bin";
			bin_fds[sync_number].open(bfname.c_str(), std::ios::in | std::ios::binary);
			if (!bin_fds[sync_number]) helper_tools::error("Cannot open file [" + bfname + "] for reading");
			//Read PED file
			std::string ped_fname = helper_tools::get_name_from_vcf(fname) + ".fam";
			std::ifstream fdp(ped_fname);
			if (!fdp.is_open()) helper_tools::error("Cannot open pedigree file [" + ped_fname + "] for reading");
			ind_names.push_back(std::vector < std::string >());
			ind_fathers.push_back(std::vector < std::string >());
			ind_mothers.push_back(std::vector < std::string >());
			while (getline(fdp, buffer)) {
				helper_tools::split(buffer, tokens);
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
		if (flagSEEK && nsamples != 0) helper_tools::error("Binary file found for a non-empty BCF file");

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

	int32_t getSamples(uint32_t file, std::map < std::string, int32_t > & samples) {
		samples.clear();
		for (uint32_t i = 0 ; i < ind_number[file] ; i ++) samples.insert(std::pair < std::string, int32_t > (ind_names[file][i], i));
		return samples.size();
	}

	//GET SAMPLE COUNTS
	int32_t getSamples() { return std::accumulate(ind_number.begin(), ind_number.end(), 0); }
	int32_t getSamples(uint32_t file) { return ind_number[file]; }

	//GET ACs/ANs/AFs
	uint32_t getAC(uint32_t file) { return AC[file]; }
	uint32_t getAN(uint32_t file) { return AN[file]; }
	uint32_t getAC() { return std::accumulate(AC.begin(), AC.end(), 0); }
	uint32_t getAN() { return std::accumulate(AN.begin(), AN.end(), 0); }
	float getAF(uint32_t file) { return AC[file]*1.0f/AN[file]; }
	float getAF() { return std::accumulate(AC.begin(), AC.end(), 0)*1.0f/std::accumulate(AN.begin(), AN.end(), 0); }



	//SET SYNCHRONIZED READER TO NEXT RECORD
	int32_t nextRecord() {

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
		for (uint32_t r = 0, firstfile = 1 ; r < sync_number ; r++) {

			//Check if reader has a record
			bool hasRecord = bcf_sr_has_line(sync_reader, r);

			//First time we see a record for this variant
			if (hasRecord) {

				//Get the record
				sync_lines[r] = bcf_sr_get_line(sync_reader, r);

				//If bi-allelic, proceed
				if (sync_lines[r]->n_allele == 2) {

					//If first time we see the record across files
					if (firstfile) {

						//Get variant information
						chr = bcf_hdr_id2name(sync_reader->readers[0].header, sync_lines[r]->rid);
						pos = sync_lines[r]->pos + 1;
						rsid = std::string(sync_lines[r]->d.id);
						ref = std::string(sync_lines[r]->d.allele[0]);
						alt = std::string(sync_lines[r]->d.allele[1]);
						firstfile = 0;
					}

					//Get AC/AN information
					int32_t rAC = bcf_get_info_int32(sync_reader->readers[r].header, sync_lines[r], "AC", &vAC, &nAC);
					int32_t rAN = bcf_get_info_int32(sync_reader->readers[r].header, sync_lines[r], "AN", &vAN, &nAN);
					if (rAC != 1) helper_tools::error("AC field is needed in file");
					if (rAN != 1) helper_tools::error("AN field is needed in file");
					AC[r] = vAC[0]; AN[r] = vAN[0];

					//Get SEEK information
					if (sync_types[r] == FILE_BINARY) {
						int32_t rsk = bcf_get_info_int32(sync_reader->readers[r].header, sync_lines[r], "SEEK", &vSK, &nSK);
						if (nSK != 4) helper_tools::error("INFO/SEEK field should contain 4 numbers");
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

					//Set "has record" flag
					sync_flags[r] = true;
				}
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

	//RETURN THE SIZE OF RECORD IN BYTES
	int32_t sizeRecord(uint32_t file) {
		return bin_size[file];
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
			int32_t rdp = bcf_get_genotypes(sync_reader->readers[file].header, sync_lines[file], buffer, &ndp);
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
	xcf_writer(std::string _hts_fname, bool _hts_genotypes, uint32_t _nthreads) {
		std::string oformat;
		hts_fname = _hts_fname;

		if (hts_fname == "-") {
			oformat = "wbu";		//Uncompressed BCF for stdout
		} else if (hts_fname.size() > 3 && hts_fname.substr(hts_fname.size()-3) == "bcf") {
			oformat = "wb";			//Compressed BCF for file
		} else if (hts_fname.size() > 5 && hts_fname.substr(hts_fname.size()-5) == "vcf.gz") {
			oformat = "wz";			//Compressed VCF for file
		} else if (hts_fname.size() > 3 && hts_fname.substr(hts_fname.size()-5) == "vcf") {
			oformat = "wv";			//Uncompressed VCF for file
		} else helper_tools::error("Filename extension of [" + hts_fname + "] not recognized");

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
			std::string bfname = helper_tools::get_name_from_vcf(hts_fname) + ".bin";
			bin_fds.open(bfname.c_str(), std::ios::out | std::ios::binary);
			if (!bin_fds) helper_tools::error("Cannot open file [" + bfname + "] for writing");
		}
	}

	//DESTRUCTOR
	~xcf_writer() {
		//close();
	}

	//Write sample IDs
	void writeHeader(std::vector < std::string > & samples, std::string contig, std::string source) {
		//
		hts_hdr = bcf_hdr_init("w");
		bcf_hdr_append(hts_hdr, std::string("##fileDate="+helper_tools::date()).c_str());
		bcf_hdr_append(hts_hdr, std::string("##source=" + source).c_str());
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
			std::string ffname = helper_tools::get_name_from_vcf(hts_fname) + ".fam";
			std::ofstream fd (ffname);
			if (!fd.is_open()) helper_tools::error("Cannot open [" + ffname + "] for writing");
			for (uint32_t i = 0 ; i < samples.size() ; i++) fd << samples[i] << "\tNA\tNA" << std::endl;
			fd.close();
		}
		if (bcf_hdr_write(hts_fd, hts_hdr) < 0) helper_tools::error("Failing to write VCF/header");
		bcf_clear1(hts_record);
	}

	//Write sample IDs
	void writeHeader(std::vector < std::string > & samples, std::vector < int > & fathers, std::vector < int > & mothers, std::string contig, std::string source) {
		//
		hts_hdr = bcf_hdr_init("w");
		bcf_hdr_append(hts_hdr, std::string("##fileDate="+helper_tools::date()).c_str());
		bcf_hdr_append(hts_hdr, std::string("##source=" + source).c_str());
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
			std::string ffname = helper_tools::get_name_from_vcf(hts_fname) + ".fam";
			std::ofstream fd (ffname);
			if (!fd.is_open()) helper_tools::error("Cannot open [" + ffname + "] for writing");
			for (uint32_t i = 0 ; i < samples.size() ; i++) {
				std::string sfather = (fathers[i]>=0)?samples[fathers[i]]:"NA";
				std::string smother = (mothers[i]>=0)?samples[mothers[i]]:"NA";
				fd << samples[i] << "\t" << sfather << "\t" << smother << std::endl;
			}
			fd.close();
		}
		if (bcf_hdr_write(hts_fd, hts_hdr) < 0) helper_tools::error("Failing to write VCF/header");
		bcf_clear1(hts_record);
	}

	//Write variant information
	void writeInfo(std::string chr, uint32_t pos, std::string ref, std::string alt, std::string rsid, uint32_t AC, uint32_t AN) {
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
		if (bcf_write1(hts_fd, hts_hdr, hts_record) < 0) helper_tools::error("Failing to write VCF/record for rare variants");
		bcf_clear1(hts_record);
	}

	void close() {
		free(vsk);
		bcf_destroy1(hts_record);
		bcf_hdr_destroy(hts_hdr);
		if (hts_close(hts_fd)) helper_tools::error("Non zero status when closing [" + hts_fname + "]");
		if (hts_fname!="-") {
			if (bcf_index_build3(hts_fname.c_str(), NULL, 14, nthreads) < 0) helper_tools::error("Fail to index file [" + hts_fname + "]");
		}

	}
};

#endif

