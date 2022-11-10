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

#include <containers/conditioning_set/conditioning_set_header.h>

using namespace std;

struct solver_callback_params {
	genotype_set * GS;
	conditioning_set * CS;
};

void * solver_callback(void * ptr) {
	solver_callback_params * P  = static_cast < solver_callback_params * >( ptr );

	int id_worker, id_job;
	pthread_mutex_lock(&P->CS->mutex_workers);
	id_worker = P->CS->i_worker ++;
	pthread_mutex_unlock(&P->CS->mutex_workers);

	for(;;) {
		pthread_mutex_lock(&P->CS->mutex_workers);
		id_job = P->CS->i_job ++;
		pthread_mutex_unlock(&P->CS->mutex_workers);
		if (id_job <= P->CS->sites_pbwt_mthreading.back()) {
			P->CS->solve(id_job, P->GS);
			pthread_mutex_lock(&P->CS->mutex_workers);
			vrb.progress("  * PBWT phasing sweep", (++P->CS->d_job)*1.0/(P->CS->sites_pbwt_mthreading.back()+1));
			pthread_mutex_unlock(&P->CS->mutex_workers);

		} else pthread_exit(NULL);
	}
}

void conditioning_set::solve(int chunk, genotype_set * GS) {

	//Allocate
	vector < int > A = vector < int > (n_hap, 0);
	vector < int > B = vector < int > (n_hap, 0);
	vector < int > C = vector < int > (n_hap, 0);
	vector < int > D = vector < int > (n_hap, 0);
	vector < int > R = vector < int > (n_hap, 0);
	vector < int > G = vector < int > (n_hap, 0);
	vector < bool > Het = vector < bool > (n_ind, 0);
	vector < bool > Mis = vector < bool > (n_ind, 0);
	vector < bool > Amb = vector < bool > (n_ind, 0);

	iota(A.begin(), A.end(), 0);
	fill(C.begin(), C.end(), 0);

	for (int l = 0 ; l < n_site ; l ++) {
		bool chnk = (sites_pbwt_mthreading[l] == chunk);
		bool buff = (sites_pbwt_mthreading[l] < chunk) && (l >= starts_pbwt_mthreading[chunk]);

		if (chnk && l) {

			//1. INIT
			double thresh = 2.5, s, s0, s1;
			unsigned int nm = 0, nh = 0;
			for (int h = 0 ; h < n_hap ; h++) G[h] = (H_opt_var.get(l, h)?1:-1);
			for (int i = 0 ; i < n_ind ; i ++) {
				Mis[i] = VAR_GET_MIS(MOD2(l), GS->vecG[i]->Variants[DIV2(l)]);
				Het[i] = VAR_GET_HET(MOD2(l), GS->vecG[i]->Variants[DIV2(l)]);
				Amb[i] = (Het[i] || Mis[i]);
				if (Amb[i]) { G[2*i+0] = 0; G[2*i+1] = 0;}
				nh+=Het[i];
			}

			//2. PHASING FIRST PASS
			while (nh && thresh > 1.0) {
				int nhOld = nh; nh = 0, nm = 0 ;
				for (int i = 0, h = 0 ; i < n_ind ; i++, h += 2) {
					if (Amb[i]) {
						if (Het[i]) {
							s = 0.0;
							if (R[h+0]>0) s += G[A[R[h+0]-1]];
							if (R[h+0]<(n_hap-1)) s += G[A[R[h+0]+1]];
							if (R[h+1]>0) s -= G[A[R[h+1]-1]];
							if (R[h+1]<(n_hap-1)) s -= G[A[R[h+1]+1]];
							if (s > thresh) { G[h+0] = 1.0; G[h+1] = -1.0; Amb[i] = false; }
							else if (s < -thresh) { G[h+0] = -1.0; G[h+1] = 1.0; Amb[i] = false; }
							else ++nh;
						}
						if (Mis[i]) {
							if (R[h+0]>0) s0 = G[A[R[h+0]-1]];
							if (R[h+0]<(n_hap-1)) s0 += G[A[R[h+0]+1]];
							if (R[h+1]>0) s1 = G[A[R[h+1]-1]];
							if (R[h+1]<(n_hap-1)) s1 += G[A[R[h+1]+1]];
							if (s0 == -2 && s1 == -2) { G[h+0] = -1.0; G[h+1] = -1.0; Amb[i] = false; }
							else if (s0 == -2 && s1 == 2) { G[h+0] = -1.0; G[h+1] = 1.0; Amb[i] = false; }
							else if (s0 == 2 && s1 == -2) { G[h+0] = 1.0; G[h+1] = -1.0; Amb[i] = false; }
							else if (s0 == 2 && s1 == 2) { G[h+0] = 1.0; G[h+1] = 1.0; Amb[i] = false; }
							else ++nm;
						}
					}
				}
				if (nh == nhOld) thresh -= 1.0 ;
			}

			//3. PHASING SECOND PASS
			if (nh || nm) {
				for (int i = 0, h = 0 ; i < n_ind ; ++i, h += 2) {
					if (Amb[i]) {
						if (Het[i]) {
							s = 0.0;
							if (R[h+0]>0) s += G[A[R[h+0]-1]] * scoreBit[l - C[R[h+0]] + 1];
							if (R[h+0]<(n_hap-1)) s += G[A[R[h+0]+1]] * scoreBit[l - C[R[h+0]+1]+1];
							if (R[h+1]>0) s -= G[A[R[h+1]-1]] * scoreBit[l - C[R[h+1]] + 1];
							if (R[h+1]<(n_hap-1)) s -= G[A[R[h+1]+1]] * scoreBit[l - C[R[h+1]+1] + 1];
							if (s > 0) { G[h+0] = 1 ; G[h+1] = -1 ; }
							else { G[h+0] = -1 ; G[h+1] = 1 ; }
						}
						if (Mis[i]) {
							if (R[h+0]>0) s0 = G[A[R[h+0]-1]] * scoreBit[l - C[R[h+0]] + 1];
							if (R[h+0]<(n_hap-1)) s0 += G[A[R[h+0]+1]] * scoreBit[l - C[R[h+0]+1]+1];
							if (R[h+1]>0) s1 = G[A[R[h+1]-1]] * scoreBit[l - C[R[h+1]] + 1];
							if (R[h+1]<(n_hap-1)) s1 += G[A[R[h+1]+1]] * scoreBit[l - C[R[h+1]+1] + 1];
							if (s0 > 0) G[h+0] = 1.0;
							else G[h+0] = -1.0;
							if (s1 > 0) G[h+1] = 1.0;
							else G[h+1] = -1.0;
						}
					}
				}
			}

			//4. UPDATE HAPS IN MATRIX
			for (int i = 0 ; i < n_ind ; i++) {
				if (Het[i] || Mis[i]) {
					H_opt_var.set(l, 2*i+0, G[2*i+0] > 0);
					H_opt_var.set(l, 2*i+1, G[2*i+1] > 0);
				}
			}
		}

		if (chnk || buff) {
			int u = 0, v = 0, p = l, q = l;
			for (int h = 0 ; h < n_hap ; h ++) {
				int alookup = A[h], dlookup = C[h];
				if (dlookup > p) p = dlookup;
				if (dlookup > q) q = dlookup;
				if (!H_opt_var.get(l, alookup)) {
					A[u] = alookup;
					C[u] = p;
					p = 0;
					u++;
				} else {
					B[v] = alookup;
					D[v] = q;
					q = 0;
					v++;
				}
			}
			std::copy(B.begin(), B.begin()+v, A.begin()+u);
			std::copy(D.begin(), D.begin()+v, C.begin()+u);

			for (int h = 0 ; h < n_hap ; h ++) R[A[h]] = h;
		}
	}
}


void conditioning_set::solve(genotype_set * GS) {
	tac.clock();
	i_worker = 0; i_job = 0, d_job = 0;

	//
	scoreBit = vector < float > (n_site, 0.0);
	for (int l = 0 ; l < n_site ; ++l) scoreBit[l] = log (l + 1.0);

	//Perform multi-threaded selection

	solver_callback_params tp;
	tp.GS = GS;
	tp.CS = this;

	vrb.progress("  * PBWT phasing sweep", 0.0f);
	if (nthread > 1) {
		for (int t = 0 ; t < nthread ; t++) pthread_create( &id_workers[t] , NULL, solver_callback, static_cast < void * > (&tp));
		for (int t = 0 ; t < nthread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int c = 0 ; c  <= sites_pbwt_mthreading.back() ; c ++) {
		solve(c, GS);
		vrb.progress("  * PBWT phasing sweep", c*1.0/(sites_pbwt_mthreading.back()+1));
	}

	//Transpose to push new haps into H hap first
	transposeHaplotypes_V2H(false, false);

	vrb.bullet("PBWT phasing sweep (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

