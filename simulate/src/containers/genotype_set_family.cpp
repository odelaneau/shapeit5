#include <containers/genotype_set.h>


genotype_set::genotype_set(haplotype_set & _H) : H(_H) {
}


genotype_set::~genotype_set() {
	IDX.clear();
	GEN.clear();
	MIS.clear();
	SID.clear();
	DIP.clear();
}

void genotype_set::fillFamily() {
	vrb.title("Assemble haplotype data");
	for (int h = 0 ; h < IDX.size() ; h ++) {
		GEN.push_back(std::vector < bool > (H.n_var, false));
		for (int v = 0 ; v < H.n_var ; v ++) GEN.back()[v] = H.HAP[v][IDX[h][0]];
		if (IDX[h].size() == 2) {
			int break_point = rng.getInt(H.n_var);
			for (int v = break_point ; v < H.n_var ; v ++) GEN.back()[v] = H.HAP[v][IDX[h][1]];
		}
	}
	vrb.bullet("Done");
}

void genotype_set::missingFamily(float frac) {
	vrb.title("Push sporadic missing");
	unsigned int n_miss = 0;
	MIS = std::vector < std::vector < bool > > (GEN.size()/2, std::vector < bool > (GEN[0].size(), false));
	for (int i = 0 ; i < MIS.size() ; i++) for (int j = 0 ; j < MIS[i].size() ; j++) {
		MIS[i][j] = (rng.getDouble() < frac);
		n_miss += MIS[i][j];
	}
	vrb.bullet("N=" + stb.str(n_miss) + " / F=" + stb.str(n_miss * 100.0 / (MIS.size() * MIS[0].size()), 3) + "%");
}

void genotype_set::includeFamily(int n_twins, int n_obs_fam, int n_obs_fam_off, int n_hid_fam, int n_hid_fam_off) {
	vrb.title("Simulate phased family samples");
	std::vector < int > O;
	for (int h = 0 ; h < H.n_hap ; h ++) O.push_back(h);
	random_shuffle(O.begin(), O.end());

	//twins
	for (int n = 0 ; n < n_twins ; n ++) {
		IDX.push_back(std::vector < int > (1, O[2*n+0]));
		IDX.push_back(std::vector < int > (1, O[2*n+1]));
		IDX.push_back(std::vector < int > (1, O[2*n+0]));
		IDX.push_back(std::vector < int > (1, O[2*n+1]));
		SID.push_back("Twin_" + stb.str(n) + "_a");
		SID.push_back("Twin_" + stb.str(n) + "_b");
	}
	if (n_twins > 0) O.erase(O.begin(), O.begin() + 2*n_twins);
	vrb.bullet(stb.str(n_twins) + " twin pairs simulated");

	//Entire families (parents + offsprings)
	for (int n = 0 ; n < n_obs_fam ; n ++) {
		IDX.push_back(std::vector < int > (1, O[4*n+0]));
		IDX.push_back(std::vector < int > (1, O[4*n+1]));
		IDX.push_back(std::vector < int > (1, O[4*n+2]));
		IDX.push_back(std::vector < int > (1, O[4*n+3]));
		SID.push_back("FamO_" + stb.str(n) + "_moth");
		SID.push_back("FamO_" + stb.str(n) + "_fath");

		for (int c = 0 ; c < n_obs_fam_off ; c ++) {
			IDX.push_back(std::vector < int > (2, 0));
			if (rng.getDouble()<0.5) {
				IDX.back()[0] = O[4*n+0];
				IDX.back()[1] = O[4*n+1];
			} else {
				IDX.back()[0] = O[4*n+1];
				IDX.back()[1] = O[4*n+0];
			}
			IDX.push_back(std::vector < int > (2, 0));
			if (rng.getDouble()<0.5) {
				IDX.back()[0] = O[4*n+2];
				IDX.back()[1] = O[4*n+3];
			} else {
				IDX.back()[0] = O[4*n+3];
				IDX.back()[1] = O[4*n+2];
			}
			SID.push_back("FamO_" + stb.str(n) + "_child_" + stb.str(c));
		}
	}
	if (n_obs_fam > 0) O.erase(O.begin(), O.begin() + 4*n_obs_fam);
	vrb.bullet(stb.str(n_obs_fam) + " families with " + stb.str(n_obs_fam_off) + " kids simulated");

	//Only kids
	for (int n = 0 ; n < n_hid_fam ; n ++) {
		for (int c = 0 ; c < n_hid_fam_off ; c ++) {
			IDX.push_back(std::vector < int > (2, 0));
			if (rng.getDouble()<0.5) {
				IDX.back()[0] = O[4*n+0];
				IDX.back()[1] = O[4*n+1];
			} else {
				IDX.back()[0] = O[4*n+1];
				IDX.back()[1] = O[4*n+0];
			}
			IDX.push_back(std::vector < int > (2, 0));
			if (rng.getDouble()<0.5) {
				IDX.back()[0] = O[4*n+2];
				IDX.back()[1] = O[4*n+3];
			} else {
				IDX.back()[0] = O[4*n+3];
				IDX.back()[1] = O[4*n+2];
			}
			SID.push_back("FamH_" + stb.str(n) + "_child_" + stb.str(c));
		}
	}
	if (n_hid_fam > 0) O.erase(O.begin(), O.begin() + 4*n_hid_fam);
	vrb.bullet(stb.str(n_hid_fam) + " groups of " + stb.str(n_hid_fam_off) + " siblings simulated");

	//Push remaining as unrelateds
	for (int o = 0 ; o < O.size() ; o += 2) {
		IDX.push_back(std::vector < int > (1, O[o+0]));
		IDX.push_back(std::vector < int > (1, O[o+1]));
		SID.push_back("Unr_" + stb.str(o/2));
	}
	vrb.bullet(stb.str(O.size()/2) + " unrelateds simulated");
	O.clear();
}
