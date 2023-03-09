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

class xcf_reader {
public:
	//HTS part
	bcf_srs_t * sync_reader;
	uint32_t sync_number;
	std::vector < bcf1_t * > sync_lines;

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
	std::vector < uint32_t > bin_flag;			//Has binary data?
	std::vector < uint32_t > bin_type;			//Type of binary record					//Integer 1 in INFO/SEEK field
	std::vector < uint64_t > bin_seek;			//Location of binary record				//Integer 2 and 3 in INFO/SEEK field
	std::vector < uint32_t > bin_size;			//Amount of binary records in bytes		//Integer 4 in INFO/SEEK field

	//
	xcf_reader(uint32_t nthreads);
	~xcf_reader();

	//
	void setRegion(std::string region);
	void addFile(std::string fname);

	//
	void readPedigree(std::string fped);

	//
	int32_t nextRecord();
	int32_t hasRecord(uint8_t file, uint8_t type);
	int32_t readRecord(uint8_t file, uint8_t type, char *);

	//
	int32_t closeFile(uint8_t ifile);
	int32_t closeFiles();
};

xcf_reader::xcf_reader(uint32_t nthreads) {
	sync_number = 0;
	sync_reader = bcf_sr_init();
	sync_reader->collapse = COLLAPSE_NONE;
	sync_reader->require_index = 1;
	if (nthreads > 1) bcf_sr_set_threads(sync_reader, nthreads);
	nsk = 0;
	rsk = 0;
	vsk = NULL;
	bcf_line = NULL;
}

xcf_reader::~xcf_reader() {
	bcf_sr_destroy(sync_reader);
}

void xcf_reader::setRegion(std::string region) {
	if (bcf_sr_set_regions(sync_reader, region.c_str(), 0) == -1)
		vrb.error("Impossible to jump to region [" + region + "]");
}

void xcf_reader::readPedigree(std::string fname) {
	std::string buffer;
	std::vector < std::string > tokens;

	//Parsing
	input_file fd(fname);
	if (fd.fail()) vrb.error("Cannot open pedigree file [" + fname + "] for reading");
	while (getline(fd_list, buffer)) {
		stb.split(buffer, tokens);
		ind_names[sync_number].push_back(tokens[0]);
		if (tokens.size() >=3 ) {
			ind_fathers[sync_number].push_back(tokens[1]);
			ind_mothers[sync_number].push_back(tokens[2]);
		} else {
			ind_fathers[sync_number].push_back("NA");
			ind_mothers[sync_number].push_back("NA");
		}
	}
	fd.close();
}

uint32_t xcf_reader::addFile(std::string fname) {
	uint32_t nbinaryfiles = 0;

	//Open BCF file and add it in sync reader
	if (!(bcf_sr_add_reader (sync_reader, fname.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf:			vrb.error("Opening [" + fname + "]: not compressed with bgzip"); break;
		case idx_load_failed: 	vrb.error("Opening [" + fname + "]: impossible to load index file"); break;
		case file_type_error: 	vrb.error("Opening [" + fname + "]: file format not supported by HTSlib"); break;
		default : 				vrb.error("Opening [" + fname + "]: unknown error"); break;
		}
	}

	//Allocation for sample informations
	ind_number.push_back(0);
	ind_names.push_back(std::vector < std::string >());
	ind_fathers.push_back(std::vector < std::string >());
	ind_mothers.push_back(std::vector < std::string >());

	//Process sample IDs
	int n_samples = bcf_hdr_nsamples(sr->readers[sync_reader].header);
	if (n_samples) {
		//Try to read directly in BCF file
		for (int i = 0 ; i < n_samples ; i ++) {
			ind_names[sync_number].push_back(std::string(sync_reader->readers[sync_number].header->samples[i]));
			ind_fathers[sync_number].push_back("NA");
			ind_mothers[sync_number].push_back("NA");
		}
	} else {
		//Open FAM / pedigree file
		std::string ped_fname = stb.get_name_from_vcf(fname) + ".fam";
		readPedigree(sync_number, fname);
	}
	ind_number[sync_number] = ind_names[sync_number].size();

	//Allocation for Binary record information
	sync_lines.push_back(NULL);
	sync_lines[sync_number] = bcf_init();
	bin_fds.push_back(std::ifstream ());
	bin_flag.push_back(0);
	bin_type.push_back(0);
	bin_seek.push_back(0);
	bin_size.push_back(0);

	//Check header for associated binary file
	int32_t id = bcf_hdr_id2int(sync_reader->readers[sync_number].header, BCF_DT_ID, "SEEK");
	if (bcf_hdr_idinfo_exists(sync_reader->readers[sync_number].header, BCF_HL_INFO, id)) {
		//There is a binary file as the INFO/SEEK is there
		std::string bfname = stb.get_name_from_vcf(fname) + ".bin";
		bin_fds[sync_number].open(bfname.c_str(), std::ios::in | std::ios::binary);
		if (!bin_fds[sync_number]) vrb.error("Cannot open file [" + bfname + "] for reading");
		bin_flag[sync_number] = 1;
	} else {
		//No binary file therefore read data from BCF
		bin_flag[sync_number] = 0;
	}

	//
	sync_number++;
}

int32_t xcf_reader::nextRecord() {
	int32_t ret = bcf_sr_next_line (sync_reader);
	if (ret) {
		for (uint32_t r = 0 ; r < sync_number ; r++) {
			if (bcf_sr_has_line(sr, r)) {
				sync_lines[sync_number] = bcf_sr_get_line(sr, 0);
				if (bin_flag[r]) {
					rsk = bcf_get_info_int32(sync_reader->readers[r].header, bcf_line, "SEEK", &vsk, &nsk);
					if (nsk != 4) vrb.error("INFO/SEEK field should contain 4 numbers");
					else {
						bin_type[r] = vsk[0];
						bin_seek[r] = vsk[1];
						bin_seek[r] *= MOD30BITS;
						bin_seek[r] += vsk[2];
						bin_size[r] = vsk[3];
						assert(vsk[3]>=0);
					}
				}
			}
		}
	}
	return ret;
}

int32_t xcf_reader::hasRecord(uint32_t file, uint8_t type) {
	return bin_flag_record[file][type];
}

int32_t xcf_reader::readRecord(uint32_t file, uint8_t type, char * buffer) {
	if (!bin_flag_record[file][type]) return -1;
	else {
		if (ifile == XCF_BCF) {
			int32_t ngt, rgt;
			bcf_line = bcf_sr_get_line(sync_reader, file);
			rgt = bcf_get_genotypes(sync_reader->readers[file].header, bcf_line, &buffer, &ngt);
			return rgt;
		} else {
			bin_fds[file][type].seekg(bin_iseek[r][f], bin_fds[file][type].beg);
			bin_fds[file][type].read(buffer, bin_nseek[r][f]);
			return bin_nseek[r][f];
		}
	}
}

#endif

