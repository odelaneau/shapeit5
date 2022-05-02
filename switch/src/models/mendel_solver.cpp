#include <models/mendel_solver.h>

mendel_solver::mendel_solver(haplotype_set & _H) : H(_H) {
	Errors = vector < vector < bool > > (H.vecSamples.size(), vector < bool > (H.n_variants, false));
	H.Phased = vector < vector < bool > > (H.vecSamples.size(), vector < bool > (H.n_variants, false));
}

mendel_solver::~mendel_solver() {
	Errors.clear();
}

void mendel_solver::set() {
	vrb.title("Assume input truth is already phased"); tac.clock();
	for (int i = 0 ; i < H.vecSamples.size() ; i++) {
		for (int l = 0 ; l < H.n_variants ; l ++) {
			int g = H.Htrue[2*i+0][l] + H.Htrue[2*i+1][l];
			H.Phased[i][l] = (g == 1);
		}
	}
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void mendel_solver::solveT(int locus, int cidx, int fidx, int midx) {
	int phased = 0, mendel = 0;
	int cg = H.Htrue[2*cidx+0][locus] + H.Htrue[2*cidx+1][locus];
	int fg = H.Htrue[2*fidx+0][locus] + H.Htrue[2*fidx+1][locus];
	int mg = H.Htrue[2*midx+0][locus] + H.Htrue[2*midx+1][locus];
	bool c0 = false, c1 = false, f0 = false, f1 = false, m0 = false, m1 = false;
	if (fg == 0 && mg == 0 && cg == 0) { f0 = 0; f1 = 0; m0 = 0; m1 = 0; c0 = 0; c1 = 0; mendel = 0; phased = 1;}
	if (fg == 0 && mg == 0 && cg == 1) { f0 = 0; f1 = 0; m0 = 0; m1 = 0; c0 = 0; c1 = 1; mendel = 1; phased = 0;}
	if (fg == 0 && mg == 0 && cg == 2) { f0 = 0; f1 = 0; m0 = 0; m1 = 0; c0 = 1; c1 = 1; mendel = 1; phased = 0;}
	if (fg == 0 && mg == 1 && cg == 0) { f0 = 0; f1 = 0; m0 = 1; m1 = 0; c0 = 0; c1 = 0; mendel = 0; phased = 1;}
	if (fg == 0 && mg == 1 && cg == 1) { f0 = 0; f1 = 0; m0 = 0; m1 = 1; c0 = 0; c1 = 1; mendel = 0; phased = 1;}
	if (fg == 0 && mg == 1 && cg == 2) { f0 = 0; f1 = 0; m0 = 0; m1 = 1; c0 = 1; c1 = 1; mendel = 1; phased = 0;}
	if (fg == 0 && mg == 2 && cg == 0) { f0 = 0; f1 = 0; m0 = 1; m1 = 1; c0 = 0; c1 = 0; mendel = 1; phased = 0;}
	if (fg == 0 && mg == 2 && cg == 1) { f0 = 0; f1 = 0; m0 = 1; m1 = 1; c0 = 0; c1 = 1; mendel = 0; phased = 1;}
	if (fg == 0 && mg == 2 && cg == 2) { f0 = 0; f1 = 0; m0 = 1; m1 = 1; c0 = 1; c1 = 1; mendel = 1; phased = 0;}
	if (fg == 1 && mg == 0 && cg == 0) { f0 = 0; f1 = 1; m0 = 0; m1 = 0; c0 = 0; c1 = 0; mendel = 0; phased = 1;}
	if (fg == 1 && mg == 0 && cg == 1) { f0 = 1; f1 = 0; m0 = 0; m1 = 0; c0 = 1; c1 = 0; mendel = 0; phased = 1;}
	if (fg == 1 && mg == 0 && cg == 2) { f0 = 1; f1 = 0; m0 = 0; m1 = 0; c0 = 1; c1 = 1; mendel = 1; phased = 0;}
	if (fg == 1 && mg == 1 && cg == 0) { f0 = 0; f1 = 1; m0 = 1; m1 = 0; c0 = 0; c1 = 0; mendel = 0; phased = 1;}
	if (fg == 1 && mg == 1 && cg == 1) { f0 = 0; f1 = 1; m0 = 0; m1 = 1; c0 = 0; c1 = 1; mendel = 0; phased = 0;}
	if (fg == 1 && mg == 1 && cg == 2) { f0 = 1; f1 = 0; m0 = 0; m1 = 1; c0 = 1; c1 = 1; mendel = 0; phased = 1;}
	if (fg == 1 && mg == 2 && cg == 0) { f0 = 0; f1 = 1; m0 = 1; m1 = 1; c0 = 0; c1 = 0; mendel = 1; phased = 0;}
	if (fg == 1 && mg == 2 && cg == 1) { f0 = 0; f1 = 1; m0 = 1; m1 = 1; c0 = 0; c1 = 1; mendel = 0; phased = 1;}
	if (fg == 1 && mg == 2 && cg == 2) { f0 = 1; f1 = 0; m0 = 1; m1 = 1; c0 = 1; c1 = 1; mendel = 0; phased = 1;}
	if (fg == 2 && mg == 0 && cg == 0) { f0 = 1; f1 = 1; m0 = 0; m1 = 0; c0 = 0; c1 = 0; mendel = 1; phased = 0;}
	if (fg == 2 && mg == 0 && cg == 1) { f0 = 1; f1 = 1; m0 = 0; m1 = 0; c0 = 1; c1 = 0; mendel = 0; phased = 1;}
	if (fg == 2 && mg == 0 && cg == 2) { f0 = 1; f1 = 1; m0 = 0; m1 = 0; c0 = 1; c1 = 1; mendel = 1; phased = 0;}
	if (fg == 2 && mg == 1 && cg == 0) { f0 = 1; f1 = 1; m0 = 0; m1 = 1; c0 = 0; c1 = 0; mendel = 1; phased = 0;}
	if (fg == 2 && mg == 1 && cg == 1) { f0 = 1; f1 = 1; m0 = 1; m1 = 0; c0 = 1; c1 = 0; mendel = 0; phased = 1;}
	if (fg == 2 && mg == 1 && cg == 2) { f0 = 1; f1 = 1; m0 = 0; m1 = 1; c0 = 1; c1 = 1; mendel = 0; phased = 1;}
	if (fg == 2 && mg == 2 && cg == 0) { f0 = 1; f1 = 1; m0 = 1; m1 = 1; c0 = 0; c1 = 0; mendel = 1; phased = 0;}
	if (fg == 2 && mg == 2 && cg == 1) { f0 = 1; f1 = 1; m0 = 1; m1 = 1; c0 = 0; c1 = 1; mendel = 1; phased = 0;}
	if (fg == 2 && mg == 2 && cg == 2) { f0 = 1; f1 = 1; m0 = 1; m1 = 1; c0 = 1; c1 = 1; mendel = 0; phased = 1;}
	H.Htrue[2*cidx+0][locus] = c0; H.Htrue[2*cidx+1][locus] = c1;
	H.Htrue[2*fidx+0][locus] = f0; H.Htrue[2*fidx+1][locus] = f1;
	H.Htrue[2*midx+0][locus] = m0; H.Htrue[2*midx+1][locus] = m1;
	H.Phased[cidx][locus] = (phased && cg==1);
	H.Phased[fidx][locus] = (phased && fg==1);
	H.Phased[midx][locus] = (phased && mg==1);
	Errors[cidx][locus] = mendel;
	Errors[fidx][locus] = mendel;
	Errors[midx][locus] = mendel;
}


void mendel_solver::solveD(int locus, int cidx, int pidx, bool father) {
	int phased = 0, mendel = 0;
	int cg = H.Htrue[2*cidx+0][locus] + H.Htrue[2*cidx+1][locus];
	int pg = H.Htrue[2*pidx+0][locus] + H.Htrue[2*pidx+1][locus];
	bool c0 = false, c1 = false, p0 = false, p1 = false;
	if (pg == 0 && cg == 0) { p0 = 0; p1 = 0; c0 = 0; c1 = 0; mendel = 0; phased = 1;}
	if (pg == 0 && cg == 1) { p0 = 0; p1 = 0; c0 = 0; c1 = 1; mendel = 0; phased = 1;}
	if (pg == 0 && cg == 2) { p0 = 0; p1 = 0; c0 = 1; c1 = 1; mendel = 1; phased = 0;}
	if (pg == 1 && cg == 0) { p0 = 0; p1 = 1; c0 = 0; c1 = 0; mendel = 0; phased = 1;}
	if (pg == 1 && cg == 1) { p0 = 0; p1 = 1; c0 = 0; c1 = 1; mendel = 0; phased = 0;}
	if (pg == 1 && cg == 2) { p0 = 1; p1 = 0; c0 = 1; c1 = 1; mendel = 0; phased = 1;}
	if (pg == 2 && cg == 0) { p0 = 1; p1 = 1; c0 = 0; c1 = 0; mendel = 1; phased = 0;}
	if (pg == 2 && cg == 1) { p0 = 1; p1 = 1; c0 = 1; c1 = 0; mendel = 0; phased = 1;}
	if (pg == 2 && cg == 2) { p0 = 1; p1 = 1; c0 = 1; c1 = 1; mendel = 0; phased = 1;}
	if (father) {
		H.Htrue[2*cidx+0][locus] = c0; H.Htrue[2*cidx+1][locus] = c1;
		H.Htrue[2*pidx+0][locus] = p0; H.Htrue[2*pidx+1][locus] = p1;
	} else {
		H.Htrue[2*cidx+0][locus] = c1; H.Htrue[2*cidx+1][locus] = c0;
		H.Htrue[2*pidx+0][locus] = p1; H.Htrue[2*pidx+1][locus] = p0;
	}
	H.Phased[cidx][locus] = (phased && cg==1);
	H.Phased[pidx][locus] = (phased && pg==1);
	Errors[cidx][locus] = mendel;
	Errors[pidx][locus] = mendel;
}

void mendel_solver::solve() {
	vrb.title("Mendel phasing using pedigrees"); tac.clock();
	for (int i = 0 ; i < H.vecSamples.size() ; i++) {
		for (int l = 0 ; l < H.n_variants ; l ++) {
			int fidx = ((H.Fathers[i] >= 0) && (!H.Missing[i][l]) && (!H.Missing[H.Fathers[i]][l]))?H.Fathers[i]:-1;
			int midx = ((H.Mothers[i] >= 0) && (!H.Missing[i][l]) && (!H.Missing[H.Mothers[i]][l]))?H.Mothers[i]:-1;
			if (fidx != -1 && midx != -1) solveT(l, i, fidx, midx);
			if (fidx == -1 && midx != -1) solveD(l, i, midx, false);
			if (fidx != -1 && midx == -1) solveD(l, i, fidx, true);
		}
	}

	unsigned int n_mendel_errors = 0;
	for (int i = 0 ; i < H.vecSamples.size() ; i++) for (int l = 0 ; l < H.n_variants ; l ++) n_mendel_errors += Errors[i][l];
	vrb.bullet("#Mendel errors = " + stb.str(n_mendel_errors));
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}


void mendel_solver::writePerSample(string fout) {
	tac.clock();
	vrb.title("Writing mendel errors per sample in [" + fout + "]");
	output_file fdo (fout);
	for (int i = 0 ; i < H.vecSamples.size() ; i++) {
		int fidx = H.Fathers[i];
		int midx = H.Mothers[i];
		int n_errors = 0, n_nmissing = 0;

		fdo << H.vecSamples[i];
		if (fidx != -1 && midx != -1) fdo << " " << H.vecSamples[fidx] << " " << H.vecSamples[midx];
		if (fidx != -1 && midx == -1) fdo << " " << H.vecSamples[fidx] << " -1";
		if (fidx == -1 && midx != -1) fdo << " -1 " << H.vecSamples[midx];
		if (fidx == -1 && midx == -1) fdo << " -1 -1";

		for (int l = 0 ; l < H.n_variants ; l ++) {
			n_errors += Errors[i][l];
			n_nmissing += !H.Missing[i][l];
		}
		fdo << " " << n_errors << " " << n_nmissing << " " << stb.str(n_errors * 100.0f / n_nmissing, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}


void mendel_solver::writePerVariant(string fout) {
	tac.clock();
	vrb.title("Writing mendel errors per variant in [" + fout + "]");
	output_file fdo (fout);
	for (int l = 0 ; l < H.n_variants ; l ++) {
		int n_errors = 0, n_nmissing = 0;
		for (int i = 0 ; i < H.vecSamples.size() ; i++) {
			n_errors += Errors[i][l];
			n_nmissing += !H.Missing[i][l];
		}
		fdo << H.RSIDs[l]  << " " << H.Positions[l] << " " << n_errors << " " << n_nmissing << " " << stb.str(n_errors * 100.0f / n_nmissing, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

