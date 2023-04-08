#include <containers/genotype_set.h>

void genotype_set::fillHaploid() {
	vrb.title("Assemble Haploid/Diploid haplotype data");
	for (int h = 0 ; h < IDX.size() ; h ++) {
		GEN.push_back(std::vector < bool > (H.n_var, false));
		for (int v = 0 ; v < H.n_var ; v ++) GEN.back()[v] = H.HAP[v][IDX[h][0]];
	}
	vrb.bullet("Done");
}

void genotype_set::missingHaploid(float fracH, float fracD) {
	vrb.title("Push sporadic missing in Haploids/Diploids");
	uint64_t n_miss_haploid = 0, n_miss_diploid = 0;
	uint64_t n_haploid = 0, n_diploid = 0;
	MIS = std::vector < std::vector < bool > > (GEN.size()/2, std::vector < bool > (GEN[0].size(), false));
	for (int i = 0 ; i < MIS.size() ; i++) {
		for (int j = 0 ; j < MIS[i].size() ; j++) {
			if (DIP[i]) {
				MIS[i][j] = (rng.getDouble() < fracD);
				n_miss_diploid += MIS[i][j];
			} else{
				MIS[i][j] = (rng.getDouble() < fracH);
				n_miss_haploid += MIS[i][j];
			}
		}
		n_haploid += !DIP[i];
		n_diploid += DIP[i];
	}
	vrb.bullet("Haploid N=" + stb.str(n_miss_haploid) + " / F=" + stb.str(n_miss_haploid * 100.0 / (n_haploid * MIS[0].size()), 3) + "%");
	vrb.bullet("Diploid N=" + stb.str(n_miss_diploid) + " / F=" + stb.str(n_miss_diploid * 100.0 / (n_diploid * MIS[0].size()), 3) + "%");
}

void genotype_set::includeHaploid(int n_haploids, int n_diploids) {
	vrb.title("Simulate phased haploid / diploid samples");
	std::vector < int > O;
	for (int h = 0 ; h < H.n_hap ; h ++) O.push_back(h);
	random_shuffle(O.begin(), O.end());

	vrb.bullet("Pool of haplotypes = " + stb.str(O.size()));

	//Diploids
	for (int n = 0 ; n < n_diploids; n ++) {
		IDX.push_back(std::vector < int > (1, O[2*n+0]));
		IDX.push_back(std::vector < int > (1, O[2*n+1]));
		SID.push_back("Dip_" + stb.str(n));
		DIP.push_back(true);
	}
	if (n_diploids > 0) O.erase(O.begin(), O.begin() + 2*n_diploids);
	vrb.bullet(stb.str(n_diploids) + " diploid samples simulated");

	vrb.bullet("Pool of haplotypes = " + stb.str(O.size()));

	//Haploids
	for (int n = 0 ; n < n_haploids; n ++) {
		IDX.push_back(std::vector < int > (1, O[n]));
		IDX.push_back(std::vector < int > (1, O[n]));
		SID.push_back("Hap_" + stb.str(n));
		DIP.push_back(false);
	}
	if (n_haploids > 0) O.erase(O.begin(), O.begin() + n_haploids);
	vrb.bullet(stb.str(n_haploids) + " haploid samples simulated");

	vrb.bullet("Pool of haplotypes = " + stb.str(O.size()));
	O.clear();
}
