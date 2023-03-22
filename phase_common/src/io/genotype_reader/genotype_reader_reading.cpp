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

#include <io/genotype_reader/genotype_reader_header.h>

#include <utils/xcf.h>
#include <containers/bitvector.h>

void genotype_reader::readGenotypes() {
	tac.clock();
	vrb.wait("  * VCF/BCF parsing");

	//File idx
	int32_t idx_file_main = 0;
	int32_t idx_file_ref = panels[1]?1:-1;
	int32_t idx_file_scaf = panels[2]?(panels[1]+panels[2]):-1;

	//Initialize VCF/BCF reader(s)
	xcf_reader XR(region, nthreads);

	//Opening file(s)
	for (int32_t f = 0 ; f < 3 ; f ++) if (panels[f]) XR.addFile(filenames[f]);

	//Main sample IDs processing
	std::vector < std::string > main_names;
	std::map < std::string, uint32_t > map_names;
	XR.getSamples(idx_file_main, main_names);
	for (int32_t i = 0 ; i < n_main_samples ; i ++) {
		G.vecG[i]->name = main_names[i];
		map_names.insert(std::pair < std::string, int32_t > (G.vecG[i]->name, i));
	}

	//Scaffold sample IDs processing
	int32_t n_scaf_samples = 0, n_with_scaffold = 0;
	std::vector < uint32_t > mappingS2M;
	if (panels[2]) {
		std::vector < std::string > scaf_names;
		n_scaf_samples = XR.getSamples(idx_file_scaf, scaf_names);
		mappingS2M = std::vector < uint32_t > (n_scaf_samples, -1);
		for (int32_t i = 0 ; i < n_scaf_samples ; i ++) {
			std::map < std::string, uint32_t > :: iterator it = map_names.find(scaf_names[i]);
			if (it != map_names.end()) { mappingS2M[i] = it->second; n_with_scaffold++; }
		}
		vrb.bullet2(stb.str(n_with_scaffold) + " samples with scaffold");
	}

	//BYTE BUFFER ALLOCATION
	int32_t * main_buffer = NULL, * ref_buffer = NULL, * scaf_buffer = NULL;
	main_buffer = (int32_t*) malloc(2 * n_main_samples * sizeof(int32_t));
	if (n_ref_samples) ref_buffer = (int32_t*) malloc(2 * n_ref_samples * sizeof(int32_t));
	if (n_scaf_samples) scaf_buffer = (int32_t*) malloc(2 * n_scaf_samples * sizeof(int32_t));

	//BIT BUFFER ALLOCATION
	bitvector main_bitvector, ref_bitvector, scaf_bitvector;
	main_bitvector.allocate(2 * n_main_samples);
	if (n_ref_samples) ref_bitvector.allocate(2 * n_ref_samples);
	if (n_scaf_samples) scaf_bitvector.allocate(2 * n_scaf_samples);

	//Parsing VCF/BCF
	uint32_t i_variant_total = 0, i_variant_kept = 0;
	while (XR.nextRecord()) {
		if (variant_mask[i_variant_total]) {

			// ====== Retrieve MAIN data ============== //
			int32_t main_type = XR.typeRecord(idx_file_main);

			// ... in BCF format
			if (main_type == RECORD_BCFVCF_GENOTYPE) {
				XR.readRecord(idx_file_main, reinterpret_cast< char** > (&main_buffer));
				for(int32_t i = 0 ; i < 2 * n_main_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(main_buffer[i+0])==1);
					bool a1 = (bcf_gt_allele(main_buffer[i+1])==1);
					bool mi = (main_buffer[i+0] == bcf_gt_missing || main_buffer[i+1] == bcf_gt_missing);
					bool he = !mi && a0 != a1;
					bool ho = !mi && a0 == a1;
					if (a0) VAR_SET_HAP0(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (a1) VAR_SET_HAP1(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (mi) VAR_SET_MIS(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (he) VAR_SET_HET(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (mi) { V.vec_pos[i_variant_kept]->cmis++; n_genotypes[3] ++; }
					else { V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1); V.vec_pos[i_variant_kept]->calt += a0+a1; n_genotypes[a0+a1] ++; }
				}
			}

			// ... in binary genotype format
			else if (main_type == RECORD_BINARY_GENOTYPE) {
				XR.readRecord(idx_file_main, reinterpret_cast< char** > (&main_bitvector.bytes));
				for(int32_t i = 0 ; i < 2 * n_main_samples ; i += 2) {
					bool a0 = main_bitvector.get(i+0);
					bool a1 = main_bitvector.get(i+1);
					bool mi = a0 && !a1;
					bool he = !mi && a0 != a1;
					bool ho = !mi && a0 == a1;
					if (a0) VAR_SET_HAP0(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (a1) VAR_SET_HAP1(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (mi) VAR_SET_MIS(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (he) VAR_SET_HET(MOD2(i_variant_kept), G.vecG[DIV2(i)]->Variants[DIV2(i_variant_kept)]);
					if (mi) { V.vec_pos[i_variant_kept]->cmis++; n_genotypes[3] ++; }
					else { V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1); V.vec_pos[i_variant_kept]->calt += a0+a1; n_genotypes[a0+a1] ++; }
				}
			}

			// ... format is unsupported
			else vrb.error("Unsupported record format [" + stb.str(main_type) + "] in [" + filenames[0] + "]");


			// ====== Retrieve REFERENCE data ============== //
			if (panels[1]) {
				int32_t ref_type = XR.typeRecord(idx_file_ref);

				// ... in BCF format
				if (ref_type == RECORD_BCFVCF_GENOTYPE) {
					XR.readRecord(idx_file_ref, reinterpret_cast< char** > (&ref_buffer));
					for(int32_t i = 0 ; i < 2 * n_ref_samples ; i += 2) {
						bool a0 = (bcf_gt_allele(ref_buffer[i+0])==1);
						bool a1 = (bcf_gt_allele(ref_buffer[i+1])==1);
						if (ref_buffer[i+0] == bcf_gt_missing || ref_buffer[i+1] == bcf_gt_missing) vrb.error("Missing genotype(s) in reference panel");
						if (!bcf_gt_is_phased(ref_buffer[i+1])) vrb.error("Unphased genotype(s) in reference panel");
						H.H_opt_hap.set(i+2*n_main_samples+0, i_variant_kept, a0);
						H.H_opt_hap.set(i+2*n_main_samples+1, i_variant_kept, a1);
						V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1);
						V.vec_pos[i_variant_kept]->calt += a0+a1;
						n_alleles[a0]++; n_alleles[a1]++;
					}
				}

				// ... in binary haplotype format
				else if (ref_type == RECORD_BINARY_HAPLOTYPE) {
					XR.readRecord(idx_file_ref, reinterpret_cast< char** > (&ref_bitvector.bytes));
					for(int32_t i = 0 ; i < 2 * n_ref_samples ; i += 2) {
						bool a0 = ref_bitvector.get(i+0);
						bool a1 = ref_bitvector.get(i+1);
						H.H_opt_hap.set(i+2*n_main_samples+0, i_variant_kept, a0);
						H.H_opt_hap.set(i+2*n_main_samples+1, i_variant_kept, a1);
						V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1);
						V.vec_pos[i_variant_kept]->calt += a0+a1;
						n_alleles[a0]++; n_alleles[a1]++;
					}
				}

				// ... in sparse haplotype format
				else if (ref_type == RECORD_SPARSE_HAPLOTYPE) {
					int32_t n_elements = XR.readRecord(idx_file_ref, reinterpret_cast< char** > (&ref_buffer)) / sizeof(int32_t);
					bool major = (XR.getAF(idx_file_ref)>0.5f);
					for(uint32_t i = 0 ; i < 2 * n_ref_samples ; i ++) H.H_opt_hap.set(i+2*n_main_samples, i_variant_kept, major);
					for(uint32_t r = 0 ; r < n_elements ; r++) H.H_opt_hap.set(ref_buffer[r]+2*n_main_samples, i_variant_kept, !major);
					if (major) {
						V.vec_pos[i_variant_kept]->cref += n_elements;
						V.vec_pos[i_variant_kept]->calt += 2 * n_ref_samples - n_elements;
						n_alleles[0] += n_elements;
						n_alleles[1] += 2 * n_ref_samples - n_elements;
					} else {
						V.vec_pos[i_variant_kept]->cref += 2 * n_ref_samples - n_elements;
						V.vec_pos[i_variant_kept]->calt += n_elements;
						n_alleles[0] += 2 * n_ref_samples - n_elements;
						n_alleles[1] += n_elements;
					}
				}

				// ... format is unsupported
				else vrb.error("Unsupported record format [" + stb.str(main_type) + "] in [" + filenames[0] + "]");
			}

			// ====== Retrieve SCAFFOLD data ============== //
			if (panels[2] && XR.hasRecord(idx_file_scaf)) {
				int32_t scaf_type = XR.typeRecord(idx_file_scaf);

				// ... in BCF format
				if (scaf_type == RECORD_BCFVCF_GENOTYPE) {
					XR.readRecord(idx_file_scaf, reinterpret_cast< char** > (&scaf_buffer));
					for(int32_t i = 0 ; i < 2 * n_scaf_samples ; i += 2) {
						int32_t ind = mappingS2M[DIV2(i)];
						if (ind >= 0) {
							bool sa0 = (bcf_gt_allele(scaf_buffer[i+0])==1);
							bool sa1 = (bcf_gt_allele(scaf_buffer[i+1])==1);
							bool sph = (bcf_gt_is_phased(scaf_buffer[i+0]) || bcf_gt_is_phased(scaf_buffer[i+1]));
							bool smi = (scaf_buffer[i+0] == bcf_gt_missing || scaf_buffer[i+1] == bcf_gt_missing);
							if ((sa0 != sa1) && !smi && sph && VAR_GET_HET(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)])) {
								VAR_SET_SCA(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
								sa0?VAR_SET_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
								sa1?VAR_SET_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
								n_genotypes[4] ++;
							}
						}
					}
				}

				// ... in binary haplotype format
				else if (scaf_type == RECORD_BINARY_HAPLOTYPE) {
					XR.readRecord(idx_file_scaf, reinterpret_cast< char** > (&scaf_bitvector.bytes));
					for(int32_t i = 0 ; i < 2 * n_scaf_samples ; i += 2) {
						int32_t ind = mappingS2M[DIV2(i)];
						if (ind >= 0) {
							bool sa0 = scaf_bitvector.get(i+0);
							bool sa1 = scaf_bitvector.get(i+1);
							if ((sa0 != sa1) && VAR_GET_HET(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)])) {
								VAR_SET_SCA(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
								sa0?VAR_SET_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
								sa1?VAR_SET_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
								n_genotypes[4] ++;
							}
						}
					}
				}

				// ... format is unsupported
				else vrb.error("Unsupported record format [" + stb.str(main_type) + "] in [" + filenames[0] + "]");

			}
			vrb.progress("  * VCF/BCF parsing", i_variant_kept * 1.0 / n_variants);
			i_variant_kept ++;
		}
		i_variant_total++;
	}
	free(main_buffer);
	if (panels[1]) free(ref_buffer);
	if (panels[2]) free(scaf_buffer);
	XR.close();

	// Report
	uint64_t n_genotypes_total = accumulate(n_genotypes.begin(), n_genotypes.begin() + 4, 0UL);
	std::string str0 = "0/0=" + stb.str(n_genotypes[0]*100.0/n_genotypes_total, 3) + "%";
	std::string str1 = "0/1=" + stb.str(n_genotypes[1]*100.0/n_genotypes_total, 3) + "%";
	std::string str2 = "1/1=" + stb.str(n_genotypes[2]*100.0/n_genotypes_total, 3) + "%";
	std::string str3 = "./.=" + stb.str(n_genotypes[3]*100.0/n_genotypes_total, 3) + "%";
	std::string str4 = "0|1=" + stb.str(n_genotypes[4]*100.0/n_genotypes[1], 3) + "%";
	std::string str5 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing done ("+str5+")");
	vrb.bullet2("Genotypes ["+str0+", "+str1+", "+str2+", "+str3+", "+str4+"]");

	if (panels[1]) {
		uint64_t n_alleles_total = n_alleles[0] + n_alleles[1];
		str0 = "0=" + stb.str(n_alleles[0]*100.0/n_alleles_total, 3) + "%";
		str1 = "1=" + stb.str(n_alleles[1]*100.0/n_alleles_total, 3) + "%";
		vrb.bullet2("Haplotypes ["+str0+", "+str1+"]");
	}
}
