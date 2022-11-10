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

#include <models/haplotype_checker.h>

using namespace std;

haplotype_checker::haplotype_checker(haplotype_set & _H, int nbins) : H(_H) {
	Errors = vector < vector < bool > > (H.IDXesti.size(), vector < bool > (H.n_variants, false));
	Checked = vector < vector < bool > > (H.IDXesti.size(), vector < bool > (H.n_variants, false));
	Calib = vector < vector < float > > (nbins, vector < float > (3, 0.0f));
}

haplotype_checker::~haplotype_checker() {
	Errors.clear();
	Checked.clear();
}

void haplotype_checker::check() {
	vrb.title("Check phasing discordances"); tac.clock();
	unsigned long int n_missed = 0, n_incorrect = 0 ;
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		for (int l_curr = 0, l_prev = -1 ; l_curr < H.n_variants ; l_curr ++) {
			bool curr_t0 = H.Htrue[2*H.IDXesti[i]+0][l_curr];
			bool curr_t1 = H.Htrue[2*H.IDXesti[i]+1][l_curr];
			bool curr_e0 = H.Hesti[2*H.IDXesti[i]+0][l_curr];
			bool curr_e1 = H.Hesti[2*H.IDXesti[i]+1][l_curr];
			bool prev_t0, prev_t1, prev_e0, prev_e1;
			bool het_check1 = (curr_e0 != curr_e1);								//Phased haplotypes are hets
			bool het_check2 = het_check1 && (curr_t0 != curr_t1);				//Validation haplotypes are hets
			bool het_check3 = het_check2 && !(H.Missing[H.IDXesti[i]][l_curr]); //Validation haplotypes are non-missing
			bool het_check4 = het_check3 && H.Phased[H.IDXesti[i]][l_curr];		//Validation haplotypes are phased
			bool het_check5 = het_check4 && H.Estimated[H.IDXesti[i]][l_curr];	//Haplotypes have been estimated there
			if (het_check5) {
				if (l_prev >= 0) {
					prev_t0 = H.Htrue[2*H.IDXesti[i]+0][l_prev];
					prev_t1 = H.Htrue[2*H.IDXesti[i]+1][l_prev];
					prev_e0 = H.Hesti[2*H.IDXesti[i]+0][l_prev];
					prev_e1 = H.Hesti[2*H.IDXesti[i]+1][l_prev];
					Errors[i][l_curr] = ((curr_t0==prev_t0) != (curr_e0==prev_e0));
					Checked[i][l_curr] = true;

					//Calibration
					if (H.Hprob[H.IDXesti[i]][l_curr]) {
						string key = stb.str(l_curr) + "_" + stb.str(H.IDXesti[i]);
						map < string, float > :: iterator itM = H.Vprob.find(key);
						if (itM != H.Vprob.end()) {
							//cout << "Found [" << key << "]" << endl;
							if (itM->second >= 0.0f && itM->second <= 1.0f) {
								int bin = itM->second * (Calib.size()-1);
								Calib[bin][0] += itM->second;
								Calib[bin][1] += Errors[i][l_curr];
								Calib[bin][2] += Checked[i][l_curr];
							} else n_incorrect ++;
						} else n_missed ++;
					}
				}
				l_prev = l_curr;

			}
		}
	}
	unsigned int n_phasing_errors = 0, n_phased_hets = 0;
	for (int i = 0 ; i < Errors.size() ; i++) for (int l = 0 ; l < Errors[i].size() ; l ++) {
		n_phasing_errors += Errors[i][l];
		n_phased_hets += Checked[i][l];
	}
	vrb.bullet("#Phasing switch error rate = " + stb.str(n_phasing_errors * 100.0f / n_phased_hets, 5));
	vrb.bullet("#missed = " + stb.str(n_missed) + " / #incorrect = " +  stb.str(n_incorrect));
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void haplotype_checker::writePerSample(string fout) {
	tac.clock();
	vrb.title("Writing phasing switch errors per sample in [" + fout + "]");
	output_file fdo (fout);
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		int n_errors = 0, n_checked = 0;
		for (int l = 0 ; l < H.n_variants ; l ++) {
			n_errors += Errors[i][l];
			n_checked += Checked[i][l];
		}
		fdo << H.vecSamples[H.IDXesti[i]] << " " << n_errors << " " << n_checked << " " << stb.str(n_errors * 100.0f / n_checked, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void haplotype_checker::writeFlipSwitchErrorPerSample(string fout) {
	tac.clock();
	vrb.title("Writing phasing flip and switch errors per sample in [" + fout + "]");
	output_file fdo (fout);
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		vector < int > HET;
		for (int l = 0 ; l < H.n_variants ; l ++) if (Checked[i][l]) HET.push_back(l);

		int n_flips = 0, n_switches = 0, n_correct = 0;
		for (int h = 2 ; h < HET.size() ; h ++) {
			int n_errors = Errors[i][HET[h-1]] + Errors[i][HET[h]];
			n_correct += (n_errors == 0);
			n_switches += (n_errors == 1);
			n_flips += (n_errors == 2);
		}
		int total = n_switches + n_flips + n_correct;
		fdo << H.vecSamples[H.IDXesti[i]] << " " << n_switches << " " << n_flips << " " << n_correct << " " << stb.str(n_switches * 100.0f / total, 2) << " " << stb.str(n_flips * 100.0f / total, 2) << " " << stb.str(n_correct * 100.0f / total, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void haplotype_checker::writePerVariant(string fout) {
	tac.clock();
	vrb.title("Writing phasing switch errors per variant in [" + fout + "]");
	output_file fdo (fout);
	for (int l = 0 ; l < H.n_variants ; l ++) {
		int n_errors = 0, n_checked = 0;
		for (int i = 0 ; i < H.IDXesti.size() ; i++) {
			n_errors += Errors[i][l];
			n_checked += Checked[i][l];
		}
		fdo << H.RSIDs[l]  << " " << H.Positions[l] << " " << n_errors << " " << n_checked << " " << stb.str(n_errors * 100.0f / n_checked, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void haplotype_checker::writePerType(string fout) {
	tac.clock();
	vrb.title("Writing phasing switch errors per variant type in [" + fout + "]");
	output_file fdo (fout);
	vector < int > b_errors = vector < int > (2, 0);
	vector < int > b_checked = vector < int > (2, 0);
	for (int l = 0 ; l < H.n_variants ; l ++) {
		bool snp = isSNP(H.REFs[l], H.ALTs[l]);
		for (int i = 0 ; i < H.IDXesti.size() ; i++) {
			b_errors[snp] += Errors[i][l];
			b_checked[snp] += Checked[i][l];
		}
	}
	for (int b = 0 ; b < b_errors.size() ; b ++) {
		if (b_checked[b] > 0)
			fdo << b << " " << b_errors[b] << " " << b_checked[b] << " " << stb.str(b_errors[b] * 100.0f / b_checked[b], 2) << endl;
		else
			fdo << b << " 0 0 0.0" << endl;
	}

	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void haplotype_checker::writePerFrequency(string fout) {
	tac.clock();
	vrb.title("Writing phasing switch errors per frequency bin in [" + fout + "]");
	int max_mac = *max_element(std::begin(H.MAC), std::end(H.MAC));
	int min_mac = *min_element(std::begin(H.MAC), std::end(H.MAC));
	int siz_mac = max_mac - min_mac + 1;
	vrb.bullet("#bins = " + stb.str(siz_mac - 1));
	output_file fdo (fout);
	vector < int > b_errors = vector < int > (siz_mac, 0);
	vector < int > b_checked = vector < int > (siz_mac, 0);
	for (int l = 0 ; l < H.n_variants ; l ++) {
		for (int i = 0 ; i < H.IDXesti.size() ; i++) {
			b_errors[H.MAC[l]-min_mac] += Errors[i][l];
			b_checked[H.MAC[l]-min_mac] += Checked[i][l];
		}
	}
	for (int b = 0 ; b < b_errors.size() ; b ++) {
		if (b_checked[b] > 0)
			fdo << b+min_mac << " " << b_errors[b] << " " << b_checked[b] << " " << stb.str(b_errors[b] * 100.0f / b_checked[b], 2) << endl;
		else
			fdo << b+min_mac << " 0 0 0.0" << endl;
	}

	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void haplotype_checker::writeBlock(string fout) {
	tac.clock();
	vrb.title("Writing correct phasing blocks per sample in [" + fout + "]");
	output_file fdo (fout);
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		fdo << H.vecSamples[H.IDXesti[i]] << " " << H.Positions[0] << endl;
		for (int l = 1 ; l < H.n_variants ; l ++) {
			if (Errors[i][l])
				fdo << H.vecSamples[H.IDXesti[i]] << " " << H.Positions[l] << endl;
		}
		fdo << H.vecSamples[H.IDXesti[i]] << " " << H.Positions.back() << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}


void haplotype_checker::writeCalibration(string fout) {
	tac.clock();

	//Write output file
	vrb.title("Writing phasing calibration in [" + fout + "]");
	output_file fdo (fout);
	for (int c = 0 ; c < Calib.size() ; c++) {
		fdo << c << " " << c * 1.0f / Calib.size() << " " << (c+1) * 1.0f / Calib.size() << " " << Calib[c][0] << " " << Calib[c][1] << " " << Calib[c][2] << endl;
	}
	fdo.close();
}
