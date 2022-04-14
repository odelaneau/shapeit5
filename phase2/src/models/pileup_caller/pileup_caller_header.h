/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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
#ifndef _MPILEUP_CALLER_H
#define _MPILEUP_CALLER_H

#include <utils/otools.h>

#include <containers/haplotype_set.h>
#include <containers/genotype_set.h>
#include <containers/variant_map.h>

inline
char getBase (int code) {
	switch (code) {
	case 1: return 'A';
	case 2: return 'C';
	case 4: return 'G';
	case 8: return 'T';
	case 15: return 'N';
	}
	return -1;
}

class het {
public:
	int idx;
	int pos;
	bool a0;
	bool a1;
	char ref;
	char alt;

	het() {
		idx = pos = -1;
		a0 = a1 = 0;
	}

	~het() {
		idx = pos = -1;
		a0 = a1 = 0;
	}

	bool empty() {
		return idx < 0;
	}

	int distance (het & h) {
		if (h.idx < 0 || idx < 0) return numeric_limits<int>::max();
		return (h.pos>pos)?(h.pos-pos):(pos-h.pos);
	}
};

class data_caller {
public:
	samFile * fp;				// File handler
	bam_hdr_t * hdr;			// File header
	hts_idx_t * idx;			// Index handler
	hts_itr_t * iter;			// NULL if a region not specified
	int min_mapQ;				// mapQ filter
	bool opened;

	data_caller (int _min_mapQ) {
		min_mapQ = _min_mapQ;
		fp = NULL; hdr = NULL; idx = NULL; iter = NULL;
		opened = false;
	}

	~data_caller() {
		close();
	}

	void open (string fbam, string fasta) {
		if (!(fp=sam_open(fbam.c_str(), "r"))) vrb.error("Cannot open file!");
		if (!(hdr=sam_hdr_read(fp))) vrb.error("Cannot parse header!");
	    if (!(idx=sam_index_load(fp, fbam.c_str()))) vrb.error("Cannot load index!");
	    if ((fasta != "") && hts_set_fai_filename(fp, fasta.c_str())) vrb.error("Cannot open fasta");
	    opened = true;
	}

	bool isOpened() { return opened; }

	void close() {
		if (hdr) bam_hdr_destroy(hdr);
		if (idx) hts_idx_destroy(idx);
		if (fp) sam_close(fp);
		if (iter) hts_itr_destroy(iter);
		fp = NULL; hdr = NULL; idx = NULL; iter = NULL;
		opened = false;
	}

	void jump(string & chr, int start, int end) {
		string region = chr + ":" + stb.str(start) + "-" + stb.str(end);
		iter = sam_itr_querys(idx, hdr, region.c_str());
		if (!iter) vrb.error("Problem jumping to region [" + region + "]");
	}

	int chr(string & _chr) {
		return bam_name2id(hdr, _chr.c_str());
	}

	int begin() {
		return iter->beg;
	}

	int end() {
		return iter->end;
	}
};

class pir {
public:
	unsigned int l_phr : 8;
	unsigned int t_phr : 8;
	unsigned int r_phr : 8;
	unsigned int l_obs : 1;
	unsigned int t_obs : 1;
	unsigned int r_obs : 1;
	unsigned int l_all : 1;
	unsigned int t_all : 1;
	unsigned int r_all : 1;
	unsigned int dummy : 2;

	pir() {
		clear();
	}

	~pir() {
		clear();
	}

	void clear() {
		l_phr = t_phr = r_phr = l_obs = t_obs = r_obs  = l_all = t_all = r_all = 0;
	}
};

class pileup_caller {
public :
	//DATA
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;

	//FASTA
	string fai_fname;
	faidx_t* fai;

	//Sequencing stats
	unsigned long int n_rhets_total;
	unsigned long int n_rhets_pired;
	unsigned long int n_bases_mismatch;
	unsigned long int n_bases_match;
	unsigned long int n_bases_lowqual;
	unsigned long int n_bases_indel;
	unsigned long int n_bases_total;
	unsigned long int n_pirs_total;
	unsigned long int n_pirs_mismatch;

	//PIRs
	map < string, pir > R;

	//PARAMETERS
	int min_mapQ;
	int min_baseQ;

	//CONSTRUCTOR/DESTRUCTOR
	pileup_caller(haplotype_set &, genotype_set &, variant_map &, int, int);
	~pileup_caller();

	//MAIN
	void loadFASTA(string);
	void queryBAM(int ind, string);

	//ROUTINES
	het getScafHetsLeft(int ind, int vr);
	het getScafHetsRight(int ind, int vr);
	bool phaseWithPIRs(int ind, int vr, het &lhet, het &thet, het &rhet);
	void parseReads(const bam_pileup1_t * v_plp, int n_plp, het & h, int side);
};

#endif
