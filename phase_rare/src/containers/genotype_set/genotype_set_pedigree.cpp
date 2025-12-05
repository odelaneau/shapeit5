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

#include <containers/genotype_set/genotype_set_header.h>

using namespace std;

void genotype_set::phasePedigrees(string fped) {
	//DATA

	mendel_error = vector < int > (n_samples, 0);
	mendel_ydone = vector < int > (n_samples, 0);
	mendel_imput = vector < int > (n_samples, 0);
	mendel_ndone = vector < int > (n_samples, 0);

	tac.clock();
	vrb.title("Trio/Duo pre-phasing");

	//Reading file
	string buffer;
	vector < string > samples_ped, fathers_ped, mothers_ped, tokens;
	input_file fd_ped(fped);
	if (fd_ped.fail()) vrb.error("Cannot open pedigree file");
	while (getline(fd_ped, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 3) vrb.error("Problem in pedigree file; each line should have at least 3 columns");
		samples_ped.push_back(tokens[0]);
		fathers_ped.push_back(tokens[1]);
		mothers_ped.push_back(tokens[2]);
	}
	fd_ped.close();
	vrb.bullet("PED file parsing [n=" + stb.str(samples_ped.size()) + "]");

	//Build map
	map < string, int > samples_str2idx;
	vector < int > fathers_idx, mothers_idx;
	for (int i = 0 ; i < n_samples ; i ++) samples_str2idx.insert(pair < string, int > (names[i], i));
	int n_trios = 0, n_duosF = 0, n_duosM = 0;
	fathers_idx = vector < int > (n_samples, -1);
	mothers_idx = vector < int > (n_samples, -1);
	for (int i = 0 ; i < samples_ped.size() ; i ++) {
		map < string, int >::iterator itK = samples_str2idx.find(samples_ped[i]);
		map < string, int >::iterator itF = samples_str2idx.find(fathers_ped[i]);
		map < string, int >::iterator itM = samples_str2idx.find(mothers_ped[i]);
		if (itK != samples_str2idx.end()) {
			fathers_idx[itK->second] = (itF != samples_str2idx.end()) ? (itF->second) : (-1) ;
			mothers_idx[itK->second] = (itM != samples_str2idx.end()) ? (itM->second) : (-1) ;
			if (fathers_idx[itK->second] >= 0 && mothers_idx[itK->second] >= 0) n_trios++;
			else if (fathers_idx[itK->second] >= 0) n_duosF ++;
			else if (mothers_idx[itK->second] >= 0) n_duosM ++;
		}
	}
	vrb.bullet("#trios = " + stb.str(n_trios) + " / #paternal_duos = " +  stb.str(n_duosF) + " / #maternal_duos = " +  stb.str(n_duosM));

	//Solving
	vector < int > mappingPlain2Sparse = vector < int >(n_samples);
	for (int vr = 0 ; vr < n_rare_variants ; vr ++) {

		//Who has data and where?
		fill(mappingPlain2Sparse.begin(), mappingPlain2Sparse.end(), -1);
		for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) mappingPlain2Sparse[GRvar_genotypes[vr][r].idx] = r;

		//Loop over samples
		for (int r = 0 ; r < GRvar_genotypes[vr].size() ; r ++) {

			unsigned int child_idx = GRvar_genotypes[vr][r].idx;
			char gen_father = -1, gen_mother = -1;

			//Is there paternal data for this sample?
			if (fathers_idx[child_idx] >= 0) {
				int father_sparse_idx = mappingPlain2Sparse[fathers_idx[child_idx]];
				gen_father = (father_sparse_idx>=0)?GRvar_genotypes[vr][father_sparse_idx].get(major_alleles[vr]):(2*major_alleles[vr]);
			}

			//Is there maternal data for this sample?
			if (mothers_idx[child_idx] >= 0) {
				int mother_sparse_idx = mappingPlain2Sparse[mothers_idx[child_idx]];
				gen_mother = (mother_sparse_idx>=0)?GRvar_genotypes[vr][mother_sparse_idx].get(major_alleles[vr]):(2*major_alleles[vr]);
			}

			//Case 1: Child genotype is missing and has data for the two parents
			if (GRvar_genotypes[vr][r].mis && gen_father >= 0 && gen_mother >= 0) {
				if (gen_father == 0 && gen_mother == 0) { GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 0; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_imput[child_idx] ++;}
				if (gen_father == 0 && gen_mother == 2) { GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 1; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_imput[child_idx] ++;}
				if (gen_father == 2 && gen_mother == 0) { GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 0; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_imput[child_idx] ++;}
				if (gen_father == 2 && gen_mother == 2) { GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 1; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_imput[child_idx] ++;}
			}

			//Case 2: Child genotype is het and has data for the two parents
			if (GRvar_genotypes[vr][r].het && gen_father >= 0 && gen_mother >= 0) {
				if (gen_father == 0 && gen_mother == 0) { mendel_error[child_idx] ++; }
				if (gen_father == 0 && gen_mother == 1) { GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 1; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_father == 0 && gen_mother == 2) { GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 1; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_father == 1 && gen_mother == 0) { GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 0; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_father == 1 && gen_mother == 1) { mendel_ndone[child_idx] ++; }
				if (gen_father == 1 && gen_mother == 2) { GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 1; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_father == 2 && gen_mother == 0) { GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 0; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_father == 2 && gen_mother == 1) { GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 0; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_father == 2 && gen_mother == 2) { mendel_error[child_idx] ++; }
			}

			//Case 3: Child genotype is het and has data only for the father
			if (GRvar_genotypes[vr][r].het && gen_father >= 0 && gen_mother < 0) {
				if (gen_father == 0) { GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 1; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_father == 1) { mendel_ndone[child_idx] ++; }
				if (gen_father == 2) { GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 0; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
			}

			//Case 4: Child genotype is het and has data only for the mother
			if (GRvar_genotypes[vr][r].het && gen_father < 0 && gen_mother >= 0) {
				if (gen_mother == 0) { GRvar_genotypes[vr][r].al0 = 1; GRvar_genotypes[vr][r].al1 = 0; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
				if (gen_mother == 1) { mendel_ndone[child_idx] ++; }
				if (gen_mother == 2) { GRvar_genotypes[vr][r].al0 = 0; GRvar_genotypes[vr][r].al1 = 1; GRvar_genotypes[vr][r].pha = 1; GRvar_genotypes[vr][r].prob = 2.0f; mendel_ydone[child_idx] ++; }
			}
		}
	}

	unsigned long n_solved = accumulate(mendel_ydone.begin(), mendel_ydone.end(), 0);
	unsigned long n_imputed = accumulate(mendel_imput.begin(), mendel_imput.end(), 0);
	unsigned long n_unsolved = accumulate(mendel_ndone.begin(), mendel_ndone.end(), 0);
	unsigned long n_errors = accumulate(mendel_error.begin(), mendel_error.end(), 0);
	unsigned long sum = n_solved + n_imputed + n_unsolved + n_errors;
	float p_solved = n_solved *100.0f / sum;
	float p_unsolved = n_unsolved *100.0f / sum;
	float p_imputed = n_imputed *100.0f / sum;
	float p_errors = n_errors *100.0f / sum;

	vrb.bullet("Solving at n=" + stb.str(sum) + " ambiguous genotypes");
	vrb.bullet2("#Hets_solved   = " + stb.str(n_solved) + " (" + stb.str(p_solved, 3) + "%)");
	vrb.bullet2("#Hets_unsolved = " + stb.str(n_unsolved) + " (" + stb.str(p_unsolved, 3) + "%)");
	vrb.bullet2("#Miss_imputed  = " + stb.str(n_imputed) + " (" + stb.str(p_imputed, 3) + "%)");
	vrb.bullet2("#Errors        = " + stb.str(n_errors) + " (" + stb.str(p_errors, 3) + "%)");

	nhets_families = n_solved;
	nmiss_families = n_imputed;
}
