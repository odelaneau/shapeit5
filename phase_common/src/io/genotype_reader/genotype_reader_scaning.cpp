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


void genotype_reader::scanGenotypes() {
	tac.clock();
	vrb.wait("  * VCF/BCF scanning");

	//File idx
	int32_t idx_file_main = 0;
	int32_t idx_file_ref = panels[1]?1:-1;
	int32_t idx_file_scaf = panels[2]?(panels[1]+panels[2]):-1;

	//Initialize VCF/BCF reader(s)
	xcf_reader XR(region, nthreads);

	//Opening file(s)
	for (int f = 0 ; f < 3 ; f ++) if (panels[f]) XR.addFile(filenames[f]);

	//Sample processing
	n_main_samples = XR.getSamples(0);
	n_ref_samples = panels[1] ? XR.getSamples(1) : 0;
	if ((n_main_samples+n_ref_samples) < 50) {
		if (n_ref_samples) vrb.error("Less than 50 samples is not enough to get reliable phasing");
		else  vrb.error("Less than 50 samples is not enough to get reliable phasing, consider using a reference panel to increase sample size");
	}

	uint32_t n_variants_noverlap = 0, n_variants_multi = 0, n_variants_notsnp = 0, n_variants_rare = 0, n_variants_nscaf = 0, n_variants_nref = 0;
	uint32_t n_variants_main_format = 0, n_variants_ref_format = 0, n_variants_scaf_format = 0;
	while (XR.nextRecord()) {

		//By defaults, we do not want the variants
		variant_mask.push_back(false);

		//See which file has a record
		bool has_main = XR.hasRecord(idx_file_main);
		bool has_ref = (panels[1] && XR.hasRecord(idx_file_ref));
		bool has_scaf = (panels[2] && XR.hasRecord(idx_file_scaf));

		//If non bi-allelic
		if ((has_main+has_ref+has_scaf)==0) { n_variants_multi++; continue; }

		//Not in reference panel
		if (panels[1] && has_main && !has_ref) { n_variants_noverlap++; continue; }

		//In reference, but not in main panel
		if (panels[1] && !has_main && has_ref) { n_variants_nref++; continue; }

		//In scaffold, but not in main panel
		if (panels[2] && !has_main && has_scaf) { n_variants_nscaf++; continue; }

		//If main panel not in proper format
		int32_t main_type = XR.typeRecord(idx_file_main);
		if ((main_type != RECORD_BCFVCF_GENOTYPE) && (main_type != RECORD_BINARY_GENOTYPE)) { n_variants_main_format ++; continue; }

		//If reference panel not in proper format
		if (panels[1]) {
			int32_t ref_type = XR.typeRecord(idx_file_ref);
			if ((ref_type != RECORD_BCFVCF_GENOTYPE) && (ref_type != RECORD_BINARY_HAPLOTYPE) && (ref_type != RECORD_SPARSE_HAPLOTYPE)) { n_variants_ref_format ++; continue; }
		}

		//If scaffold panel not in proper format
		if (panels[2]) {
			int32_t scaf_type = XR.typeRecord(idx_file_scaf);
			if ((scaf_type != RECORD_BCFVCF_GENOTYPE) && (scaf_type != RECORD_BINARY_HAPLOTYPE)) { n_variants_scaf_format ++; continue; }
		}

		//Keep SNPs only
		if (filter_snp_only) {
			bool bref = (XR.ref == "A") || (XR.ref == "T") || (XR.ref == "G") || (XR.ref == "C");
			bool balt = (XR.alt == "A") || (XR.alt == "T") || (XR.alt == "G") || (XR.alt == "C");
			n_variants_notsnp += (!bref || !balt);
			if (!bref || !balt) continue;
		}

		//Keep common only
		if (filter_min_maf > 0) {
			float maf = std::min(1.0f - XR.getAF(0), XR.getAF(0));
			n_variants_rare += (maf < filter_min_maf);
			if (maf < filter_min_maf) continue;
		}

		//Push variant information
		V.push(new variant (XR.chr, XR.pos, XR.rsid, XR.ref, XR.alt, V.size()));

		//Flag it!
		variant_mask.back() = true;
		n_variants++;
	}
	XR.close();

	vrb.bullet("VCF/BCF scanning done (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("Samples [#target=" + stb.str(n_main_samples) + " / #reference=" + stb.str(n_ref_samples) + " / #sites=" + stb.str(n_variants) + " / region=" + region + "]");
	if (n_variants_noverlap) vrb.bullet2(stb.str(n_variants_noverlap) + " sites removed in main panel [not in reference panel]");
	if (n_variants_multi) vrb.bullet2(stb.str(n_variants_multi) + " sites removed in main panel [multi-allelic]");
	if (n_variants_notsnp) vrb.bullet2(stb.str(n_variants_notsnp) + " sites removed in main panel [not SNPs]");
	if (n_variants_rare) vrb.bullet2(stb.str(n_variants_rare) + " sites removed in main panel [below MAF threshold]");
	if (n_variants_nref) vrb.bullet2(stb.str(n_variants_nref) + " sites removed in reference panel [not in main panel]");
	if (n_variants_nscaf) vrb.bullet2(stb.str(n_variants_nscaf) + " sites removed in scaffold panel [not in main panel]");
	if (n_variants_main_format) vrb.bullet2(stb.str(n_variants_main_format) + " sites removed [record in main panel not in a supported format]");
	if (n_variants_ref_format) vrb.bullet2(stb.str(n_variants_ref_format) + " sites removed [record in reference panel not in a supported format]");
	if (n_variants_scaf_format) vrb.bullet2(stb.str(n_variants_scaf_format) + " sites removed [record in scaffold panel not in a supported format]");
	if (n_variants == 0) vrb.error("No variants to be phased!");
}

