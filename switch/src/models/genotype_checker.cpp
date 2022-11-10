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

#include <models/genotype_checker.h>

using namespace std;

genotype_checker::genotype_checker(haplotype_set & _H) : H(_H) {
	Errors = vector < vector < bool > > (H.IDXesti.size(), vector < bool > (H.n_variants, false));
}

genotype_checker::~genotype_checker() {
	Errors.clear();
}

void genotype_checker::check() {
	tac.clock();
	vrb.title("Check genotyping discordances");
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		for (int l = 0 ; l < H.n_variants ; l ++) {
			Errors[i][l] = (!H.Missing[H.IDXesti[i]][l]) && ((H.Htrue[2*H.IDXesti[i]+0][l] + H.Htrue[2*H.IDXesti[i]+1][l]) != (H.Hesti[2*H.IDXesti[i]+0][l] + H.Hesti[2*H.IDXesti[i]+1][l]));
		}
	}

	unsigned int n_genotyping_errors = 0;
	for (int i = 0 ; i < Errors.size() ; i++) for (int l = 0 ; l < Errors[i].size() ; l ++) n_genotyping_errors += Errors[i][l];
	vrb.bullet("#Genotyping errors = " + stb.str(n_genotyping_errors));
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void genotype_checker::writePerSample(string fout) {
	tac.clock();
	vrb.title("Writing genotyping discordances per sample in [" + fout + "]");
	output_file fdo (fout);
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		int n_errors = 0, n_nmissing = 0;
		for (int l = 0 ; l < H.n_variants ; l ++) {
			n_errors += Errors[i][l];
			n_nmissing += !H.Missing[H.IDXesti[i]][l];
		}
		fdo << H.vecSamples[H.IDXesti[i]] << " " << n_errors << " " << n_nmissing << " " << stb.str(n_errors * 100.0f / n_nmissing, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}


void genotype_checker::writePerVariant(string fout) {
	tac.clock();
	vrb.title("Writing genotyping discordances per variant in [" + fout + "]");
	output_file fdo (fout);
	for (int l = 0 ; l < H.n_variants ; l ++) {
		int n_errors = 0, n_nmissing = 0;
		for (int i = 0 ; i < H.IDXesti.size() ; i++) {
			n_errors += Errors[i][l];
			n_nmissing += !H.Missing[H.IDXesti[i]][l];
		}
		fdo << H.RSIDs[l]  << " " << H.Positions[l] << " " << n_errors << " " << n_nmissing << " " << stb.str(n_errors * 100.0f / n_nmissing, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}
