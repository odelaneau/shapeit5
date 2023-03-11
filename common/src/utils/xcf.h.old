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

#include <otools.h>

#define DATA_VOID	0
#define DATA_BCF	1
#define DATA_BINARY	2

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
	std::vector < uint32_t > sync_types;	//Type of data? VOID [BCF n=0], BCF [BCF n>0] or BINARY [BIN n>0]
	std::vector < uint32_t > sync_flags;	//Has data for the current variant?

	//Variant information
	std::string chr;
	uint32_t pos;
	std::string ref;
	std::string alt;
	std::string rsid;

	//INFO field
	int32_t nsk;
	int32_t rsk;
	int32_t * vsk;

	//Pedigree file
	std::vector < uint32_t> ind_number;
	std::vector < std::vector < std::string > > ind_names;
	std::vector < std::vector < std::string > > ind_fathers;
	std::vector < std::vector < std::string > > ind_mothers;

	//Binary files [files x types]
	std::vector < std::ifstream > bin_fds;		//File Descriptors
	std::vector < uint32_t > bin_type;			//Type of Binary record					//Integer 1 in INFO/SEEK field
	std::vector < uint64_t > bin_seek;			//Location of Binary record				//Integer 2 and 3 in INFO/SEEK field
	std::vector < uint32_t > bin_size;			//Amount of Binary records in bytes		//Integer 4 in INFO/SEEK field

	//CONSTRUCTOR
	xcf_reader::xcf_reader(uint32_t nthreads) {
		sync_number = 0;
		sync_reader = bcf_sr_init();
		sync_reader->collapse = COLLAPSE_NONE;
		sync_reader->require_index = 1;
		if (nthreads > 1) bcf_sr_set_threads(sync_reader, nthreads);
		nsk = 0;
		rsk = 0;
		vsk = NULL;
	}

	//DESTRUCTOR
	xcf_reader::~xcf_reader() {
		bcf_sr_destroy(sync_reader);
	}

	//SET REGION TO PARSE
	void xcf_reader::setRegion(std::string region) {
		if (bcf_sr_set_regions(sync_reader, region.c_str(), 0) == -1)
			vrb.error("Impossible to jump to region [" + region + "]");
	}

	//ADD A NEW FILE IN THE SYNCHRONIZED READER
	uint32_t xcf_reader::addFile(std::string fname) {
		std::string buffer;
		std::vector < std::string > tokens;

		//Open BCF file and add it in the synchronized reader
		if (!(bcf_sr_add_reader (sync_reader, fname.c_str()))) {
			switch (sr->errnum) {
			case not_bgzf:			vrb.error("Opening [" + fname + "]: not compressed with bgzip"); break;
			case idx_load_failed: 	vrb.error("Opening [" + fname + "]: impossible to load index file"); break;
			case file_type_error: 	vrb.error("Opening [" + fname + "]: file format not supported by HTSlib"); break;
			default : 				vrb.error("Opening [" + fname + "]: unknown error"); break;
			}
		}

		//Allocation for Binary record information
		sync_lines.push_back(bcf_init());
		sync_flags.push_back(0);
		bin_fds.push_back(std::ifstream());
		bin_type.push_back(0);
		bin_seek.push_back(0);
		bin_size.push_back(0);

		//Check header for associated binary file
		int32_t flagSEEK = bcf_hdr_idinfo_exists(sync_reader->readers[sync_number].header, BCF_HL_INFO, bcf_hdr_id2int(sync_reader->readers[sync_number].header, BCF_DT_ID, "SEEK"));

		//Check header for number of samples in BCF
		uint32_t nsamples = bcf_hdr_nsamples(sr->readers[sync_reader].header);

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
			while (getline(fdp, buffer)) {
				stb.split(buffer, tokens);
				ind_names[sync_number].push_back(tokens[0]);
				if (tokens.size() >=3 ) { ind_fathers[sync_number].push_back(tokens[1]); ind_mothers[sync_number].push_back(tokens[2]); }
				else { ind_fathers[sync_number].push_back("NA"); ind_mothers[sync_number].push_back("NA"); }
			}
			ind_number.push_back(ind_names[sync_number].size());
			sync_types.push_back(DATA_BINARY);
			fdp.close();
		}

		/************************************************************************************/
		/*   CASE2: There is NOT a binary file and data in the BCF 							*/
		/************************************************************************************/
		if (!flagSEEK && nsamples != 0) {
			ind_names.push_back(std::vector < std::string >());
			for (int i = 0 ; i < n_samples ; i ++)
				ind_names[sync_number].push_back(std::string(sync_reader->readers[sync_number].header->samples[i]));
			ind_number.push_back(ind_names[sync_number].size());
			ind_fathers.push_back(std::vector < std::string >("NA", ind_number[sync_number]));
			ind_mothers.push_back(std::vector < std::string >("NA", ind_number[sync_number]));
			sync_types.push_back(DATA_BCF);
		}

		/************************************************************************************/
		/*   CASE3: There is NOT a binary file and No data in the BCF 						*/
		/************************************************************************************/
		if (!flagSEEK && nsamples == 0) {
			ind_number.push_back(0);
			sync_types.push_back(DATA_VOID);
		}

		//Increment number of readers
		sync_number++;
	}

	//SET SYNCHRONIZED READER TO NEXT RECORD
	int32_t xcf_reader::nextRecord() {
		int32_t ret = bcf_sr_next_line (sync_reader);
		if (ret) {
			for (uint32_t r = 0, i = 0 ; r < sync_number ; r++) {
				sync_flags[r] = bcf_sr_has_line(sr, r);
				if (sync_flags[r]) {
					sync_lines[sync_number] = bcf_sr_get_line(sync_reader, 0);
					//Get variant information
					if (i == 0) {
						chr = bcf_hdr_id2name(sr->readers[0].header, sync_lines[sync_number]->rid);
						pos = sync_lines[sync_number]->pos + 1;
						rsid = std::string(sync_lines[sync_number]->d.id);
						ref = std::string(sync_lines[sync_number]->d.allele[0]);
						alt = std::string(sync_lines[sync_number]->d.allele[1]);
						i = 1;
					}
					//Get SEEK INFO field
					if (sync_types[r] == DATA_BINARY) {
						rsk = bcf_get_info_int32(sync_reader->readers[r].header, bcf_line, "SEEK", &vsk, &nsk);
						if (nsk != 4) vrb.error("INFO/SEEK field should contain 4 numbers");
						else {
							bin_type[r] = vsk[0];
							bin_seek[r] = vsk[1];
							bin_seek[r] *= MOD30BITS;
							bin_seek[r] += vsk[2];
							bin_size[r] = vsk[3];
						}
					}
				}
			}
		}
		return ret;
	}

	//CHECK IF FILE HAS A RECORD THERE
	int32_t xcf_reader::hasRecord(uint32_t file) {
		return (!sync_flags[file]);
	}

	//CHECK THE TYPE OF RECORD AVAILABLE FOR FILE
	// -1: No record available
	//  0: There is a record but no sample data available
	//  1: There is a record with sample data in BCF
	// >2: There is a record with sample data in Binary file / return the type of binary record
	int32_t xcf_reader::typeRecord(uint32_t file) {
		if (!sync_flags[file]) return -1;
		else if (sync_types[file] == DATA_VOID) return 0;
		else if (sync_types[file] == DATA_BCF) return 1;
		else return bin_type[r];
	}

	//READ DATA OF THE AVAILABLE RECORD
	// -1: No record available
	//  0: No sample data available
	// >0: Amount of data read in bytes
	int32_t xcf_reader::readRecord(uint32_t file, char ** buffer, uint32_t * sbuffer) {
		//No data in this file for the current record
		if (!sync_flags[file]) return -1;
		//No data since file is empty
		if (sync_types[file] == DATA_VOID) return 0;
		//Data is in BCF file directly
		if (sync_types[file] == DATA_BCF) {
			//Reallocate buffer manually so that HTSlib does not need to do it
			if (*sbuffer < (ind_number[file] * 2 * sizeof(int32_t))) {
				*buffer = (char*)realloc(*buffer, ind_number[file] * 2 * sizeof(int32_t));
				*sbuffer = ind_number[file] * 2 * sizeof(int32_t);
			}
			bcf_get_genotypes(sync_reader->readers[file].header, sync_lines[sync_number], buffer, sbuffer);
			return (ind_number[file] * 2 * sizeof(int32_t));
		}
		//Data is in Binary format
		if (bin_size[r] > *sbuffer) {
			*buffer = (char*)realloc(*buffer, bin_size[r]);
			*sbuffer = bin_size[r];
		}
		//Seek to position in Binary file
		bin_fds[file].seekg(bin_seek[file], bin_fds[file].beg);
		//Read data in Binary file
		bin_fds[file].read(*buffer, bin_size[file]);
		//Return amount of data read in file
		return bin_size[file];
	}

	int32_t xcf_reader::closeFiles() {
		free(vsk);
		for (uint32_t r = 0 ; r < sync_number ; r++) if (bin_flag[r]) bin_fds.close();
		bcf_sr_destroy(sr);
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
	htsFile * hts_fd;
	bcf1_t * hts_record;
	uint32_t hts_type;							//Type of data? VOID [BCF n=0], BCF [BCF n>0] or BINARY [BIN n>0]

	//Variant information
	std::string chr;
	uint32_t pos;
	std::string ref;
	std::string alt;
	std::string rsid;

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
	xcf_writer::xcf_writer(uint32_t nthreads) {
		sync_number = 0;
		sync_reader = bcf_sr_init();
		sync_reader->collapse = COLLAPSE_NONE;
		sync_reader->require_index = 1;
		if (nthreads > 1) bcf_sr_set_threads(sync_reader, nthreads);
		nsk = 0;
		rsk = 0;
		vsk = NULL;
	}

	//DESTRUCTOR
	xcf_reader::~xcf_reader() {
		bcf_sr_destroy(sync_reader);
	}

	//SET REGION TO PARSE
	void xcf_reader::setRegion(std::string region) {
		if (bcf_sr_set_regions(sync_reader, region.c_str(), 0) == -1)
			vrb.error("Impossible to jump to region [" + region + "]");
	}

	//ADD A NEW FILE IN THE SYNCHRONIZED READER
	uint32_t xcf_reader::addFile(std::string fname) {
		std::string buffer;
		std::vector < std::string > tokens;

		//Open BCF file and add it in the synchronized reader
		if (!(bcf_sr_add_reader (sync_reader, fname.c_str()))) {
			switch (sr->errnum) {
			case not_bgzf:			vrb.error("Opening [" + fname + "]: not compressed with bgzip"); break;
			case idx_load_failed: 	vrb.error("Opening [" + fname + "]: impossible to load index file"); break;
			case file_type_error: 	vrb.error("Opening [" + fname + "]: file format not supported by HTSlib"); break;
			default : 				vrb.error("Opening [" + fname + "]: unknown error"); break;
			}
		}

		//Allocation for Binary record information
		sync_lines.push_back(bcf_init());
		sync_flags.push_back(0);
		bin_fds.push_back(std::ifstream());
		bin_type.push_back(0);
		bin_seek.push_back(0);
		bin_size.push_back(0);

		//Check header for associated binary file
		int32_t flagSEEK = bcf_hdr_idinfo_exists(sync_reader->readers[sync_number].header, BCF_HL_INFO, bcf_hdr_id2int(sync_reader->readers[sync_number].header, BCF_DT_ID, "SEEK"));

		//Check header for number of samples in BCF
		uint32_t nsamples = bcf_hdr_nsamples(sr->readers[sync_reader].header);

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
			while (getline(fdp, buffer)) {
				stb.split(buffer, tokens);
				ind_names[sync_number].push_back(tokens[0]);
				if (tokens.size() >=3 ) { ind_fathers[sync_number].push_back(tokens[1]); ind_mothers[sync_number].push_back(tokens[2]); }
				else { ind_fathers[sync_number].push_back("NA"); ind_mothers[sync_number].push_back("NA"); }
			}
			ind_number.push_back(ind_names[sync_number].size());
			sync_types.push_back(DATA_BINARY);
			fdp.close();
		}

		/************************************************************************************/
		/*   CASE2: There is NOT a binary file and data in the BCF 							*/
		/************************************************************************************/
		if (!flagSEEK && nsamples != 0) {
			ind_names.push_back(std::vector < std::string >());
			for (int i = 0 ; i < n_samples ; i ++)
				ind_names[sync_number].push_back(std::string(sync_reader->readers[sync_number].header->samples[i]));
			ind_number.push_back(ind_names[sync_number].size());
			ind_fathers.push_back(std::vector < std::string >("NA", ind_number[sync_number]));
			ind_mothers.push_back(std::vector < std::string >("NA", ind_number[sync_number]));
			sync_types.push_back(DATA_BCF);
		}

		/************************************************************************************/
		/*   CASE3: There is NOT a binary file and No data in the BCF 						*/
		/************************************************************************************/
		if (!flagSEEK && nsamples == 0) {
			ind_number.push_back(0);
			sync_types.push_back(DATA_VOID);
		}

		//Increment number of readers
		sync_number++;
	}

	//SET SYNCHRONIZED READER TO NEXT RECORD
	int32_t xcf_reader::nextRecord() {
		int32_t ret = bcf_sr_next_line (sync_reader);
		if (ret) {
			for (uint32_t r = 0, i = 0 ; r < sync_number ; r++) {
				sync_flags[r] = bcf_sr_has_line(sr, r);
				if (sync_flags[r]) {
					sync_lines[sync_number] = bcf_sr_get_line(sync_reader, 0);
					//Get variant information
					if (i == 0) {
						chr = bcf_hdr_id2name(sr->readers[0].header, sync_lines[sync_number]->rid);
						pos = sync_lines[sync_number]->pos + 1;
						rsid = std::string(sync_lines[sync_number]->d.id);
						ref = std::string(sync_lines[sync_number]->d.allele[0]);
						alt = std::string(sync_lines[sync_number]->d.allele[1]);
						i = 1;
					}
					//Get SEEK INFO field
					if (sync_types[r] == DATA_BINARY) {
						rsk = bcf_get_info_int32(sync_reader->readers[r].header, bcf_line, "SEEK", &vsk, &nsk);
						if (nsk != 4) vrb.error("INFO/SEEK field should contain 4 numbers");
						else {
							bin_type[r] = vsk[0];
							bin_seek[r] = vsk[1];
							bin_seek[r] *= MOD30BITS;
							bin_seek[r] += vsk[2];
							bin_size[r] = vsk[3];
						}
					}
				}
			}
		}
		return ret;
	}

	//CHECK IF FILE HAS A RECORD THERE
	int32_t xcf_reader::hasRecord(uint32_t file) {
		return (!sync_flags[file]);
	}

	//CHECK THE TYPE OF RECORD AVAILABLE FOR FILE
	// -1: No record available
	//  0: There is a record but no sample data available
	//  1: There is a record with sample data in BCF
	// >2: There is a record with sample data in Binary file / return the type of binary record
	int32_t xcf_reader::typeRecord(uint32_t file) {
		if (!sync_flags[file]) return -1;
		else if (sync_types[file] == DATA_VOID) return 0;
		else if (sync_types[file] == DATA_BCF) return 1;
		else return bin_type[r];
	}

	//READ DATA OF THE AVAILABLE RECORD
	// -1: No record available
	//  0: No sample data available
	// >0: Amount of data read in bytes
	int32_t xcf_reader::readRecord(uint32_t file, char ** buffer, uint32_t * sbuffer) {
		//No data in this file for the current record
		if (!sync_flags[file]) return -1;
		//No data since file is empty
		if (sync_types[file] == DATA_VOID) return 0;
		//Data is in BCF file directly
		if (sync_types[file] == DATA_BCF) {
			//Reallocate buffer manually so that HTSlib does not need to do it
			if (*sbuffer < (ind_number[file] * 2 * sizeof(int32_t))) {
				*buffer = (char*)realloc(*buffer, ind_number[file] * 2 * sizeof(int32_t));
				*sbuffer = ind_number[file] * 2 * sizeof(int32_t);
			}
			bcf_get_genotypes(sync_reader->readers[file].header, sync_lines[sync_number], buffer, sbuffer);
			return (ind_number[file] * 2 * sizeof(int32_t));
		}
		//Data is in Binary format
		if (bin_size[r] > *sbuffer) {
			*buffer = (char*)realloc(*buffer, bin_size[r]);
			*sbuffer = bin_size[r];
		}
		//Seek to position in Binary file
		bin_fds[file].seekg(bin_seek[file], bin_fds[file].beg);
		//Read data in Binary file
		bin_fds[file].read(*buffer, bin_size[file]);
		//Return amount of data read in file
		return bin_size[file];
	}

	int32_t xcf_reader::closeFiles() {
		free(vsk);
		for (uint32_t r = 0 ; r < sync_number ; r++) if (bin_flag[r]) bin_fds.close();
		bcf_sr_destroy(sr);
	}
};

#endif

