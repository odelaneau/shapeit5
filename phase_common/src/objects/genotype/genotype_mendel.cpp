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

#include <objects/genotype/genotype_header.h>

using namespace std;

//counts[0] : # observed mendel errors
//counts[1] : # possible mendel errors
//counts[2] : # hets being scaffolded
//counts[3] : # hets not being scaffolded
void genotype::scaffoldTrio(genotype * gfather, genotype * gmother, vector < unsigned int > & counts) {
	for (int v = 0 ; v < n_variants ; v ++) {
		if (VAR_GET_HET(MOD2(v), Variants[DIV2(v)])) {
			bool father_is_hom = VAR_GET_HOM(MOD2(v), gfather->Variants[DIV2(v)]);
			bool mother_is_hom = VAR_GET_HOM(MOD2(v), gmother->Variants[DIV2(v)]);
			if (father_is_hom && mother_is_hom) {
				bool fath0 = VAR_GET_HAP0(MOD2(v), gfather->Variants[DIV2(v)]);
				bool moth0 = VAR_GET_HAP0(MOD2(v), gmother->Variants[DIV2(v)]);
				if (fath0 != moth0) {
					VAR_SET_SCA(MOD2(v), Variants[DIV2(v)]); counts[2]++;
					fath0?VAR_SET_HAP0(MOD2(v), Variants[DIV2(v)]):VAR_CLR_HAP0(MOD2(v), Variants[DIV2(v)]);
					fath0?VAR_CLR_HAP1(MOD2(v), Variants[DIV2(v)]):VAR_SET_HAP1(MOD2(v), Variants[DIV2(v)]);
				} else counts[0]++;
			} else if (father_is_hom) {
				bool fath0 = VAR_GET_HAP0(MOD2(v), gfather->Variants[DIV2(v)]);
				VAR_SET_SCA(MOD2(v), Variants[DIV2(v)]); counts[2]++;
				fath0?VAR_SET_HAP0(MOD2(v), Variants[DIV2(v)]):VAR_CLR_HAP0(MOD2(v), Variants[DIV2(v)]);
				fath0?VAR_CLR_HAP1(MOD2(v), Variants[DIV2(v)]):VAR_SET_HAP1(MOD2(v), Variants[DIV2(v)]);
			} else if (mother_is_hom) {
				bool moth0 = VAR_GET_HAP0(MOD2(v), gmother->Variants[DIV2(v)]);
				VAR_SET_SCA(MOD2(v), Variants[DIV2(v)]); counts[2]++;
				moth0?VAR_CLR_HAP0(MOD2(v), Variants[DIV2(v)]):VAR_SET_HAP0(MOD2(v), Variants[DIV2(v)]);
				moth0?VAR_SET_HAP1(MOD2(v), Variants[DIV2(v)]):VAR_CLR_HAP1(MOD2(v), Variants[DIV2(v)]);
			} else counts[3]++;
			counts[1] ++;
		} else if (VAR_GET_HOM(MOD2(v), Variants[DIV2(v)])) {
			bool father_is_hom = VAR_GET_HOM(MOD2(v), gfather->Variants[DIV2(v)]);
			bool mother_is_hom = VAR_GET_HOM(MOD2(v), gmother->Variants[DIV2(v)]);
			bool fath0 = VAR_GET_HAP0(MOD2(v), gfather->Variants[DIV2(v)]);
			bool moth0 = VAR_GET_HAP0(MOD2(v), gmother->Variants[DIV2(v)]);
			bool child0 = VAR_GET_HAP0(MOD2(v), Variants[DIV2(v)]);
			if (father_is_hom && fath0 != child0) counts[0]++;
			if (mother_is_hom && moth0 != child0) counts[0]++;
			counts[1] ++;
		} else if (VAR_GET_MIS(MOD2(v), Variants[DIV2(v)])) {
			bool father_is_hom = VAR_GET_HOM(MOD2(v), gfather->Variants[DIV2(v)]);
			bool mother_is_hom = VAR_GET_HOM(MOD2(v), gmother->Variants[DIV2(v)]);
			bool fath0 = VAR_GET_HAP0(MOD2(v), gfather->Variants[DIV2(v)]);
			bool moth0 = VAR_GET_HAP0(MOD2(v), gmother->Variants[DIV2(v)]);
			if (fath0 != moth0) VAR_SET_SCA(MOD2(v), Variants[DIV2(v)]);
			else VAR_SET_HOM(MOD2(v), Variants[DIV2(v)]);
			fath0?VAR_SET_HAP0(MOD2(v), Variants[DIV2(v)]):VAR_CLR_HAP0(MOD2(v), Variants[DIV2(v)]);
			moth0?VAR_SET_HAP1(MOD2(v), Variants[DIV2(v)]):VAR_CLR_HAP1(MOD2(v), Variants[DIV2(v)]);
		}
	}
}

void genotype::scaffoldDuoFather(genotype * gfather, vector < unsigned int > & counts) {
	for (int v = 0 ; v < n_variants ; v ++) {
		if (VAR_GET_HET(MOD2(v), Variants[DIV2(v)])) {
			bool father_is_hom = VAR_GET_HOM(MOD2(v), gfather->Variants[DIV2(v)]);
			if (father_is_hom) {
				bool fath0 = VAR_GET_HAP0(MOD2(v), gfather->Variants[DIV2(v)]);
				VAR_SET_SCA(MOD2(v), Variants[DIV2(v)]); counts[2]++;
				fath0?VAR_SET_HAP0(MOD2(v), Variants[DIV2(v)]):VAR_CLR_HAP0(MOD2(v), Variants[DIV2(v)]);
				fath0?VAR_CLR_HAP1(MOD2(v), Variants[DIV2(v)]):VAR_SET_HAP1(MOD2(v), Variants[DIV2(v)]);
			} else counts[3]++;
			counts[1] ++;
		} else if (VAR_GET_HOM(MOD2(v), Variants[DIV2(v)])) {
			bool father_is_hom = VAR_GET_HOM(MOD2(v), gfather->Variants[DIV2(v)]);
			bool fath0 = VAR_GET_HAP0(MOD2(v), gfather->Variants[DIV2(v)]);
			bool child0 = VAR_GET_HAP0(MOD2(v), Variants[DIV2(v)]);
			if (father_is_hom && fath0 != child0) counts[0]++;
			counts[1] ++;
		}
	}
}

void genotype::scaffoldDuoMother(genotype * gmother, vector < unsigned int > & counts) {
	for (int v = 0 ; v < n_variants ; v ++) {
		if (VAR_GET_HET(MOD2(v), Variants[DIV2(v)])) {
			bool mother_is_hom = VAR_GET_HOM(MOD2(v), gmother->Variants[DIV2(v)]);
			if (mother_is_hom) {
				bool moth0 = VAR_GET_HAP0(MOD2(v), gmother->Variants[DIV2(v)]);
				VAR_SET_SCA(MOD2(v), Variants[DIV2(v)]); counts[2]++;
				moth0?VAR_CLR_HAP0(MOD2(v), Variants[DIV2(v)]):VAR_SET_HAP0(MOD2(v), Variants[DIV2(v)]);
				moth0?VAR_SET_HAP1(MOD2(v), Variants[DIV2(v)]):VAR_CLR_HAP1(MOD2(v), Variants[DIV2(v)]);
			} else counts[3]++;
			counts[1] ++;
		} else if (VAR_GET_HOM(MOD2(v), Variants[DIV2(v)])) {
			bool mother_is_hom = VAR_GET_HOM(MOD2(v), gmother->Variants[DIV2(v)]);
			bool moth0 = VAR_GET_HAP0(MOD2(v), gmother->Variants[DIV2(v)]);
			bool child0 = VAR_GET_HAP0(MOD2(v), Variants[DIV2(v)]);
			if (mother_is_hom && moth0 != child0) counts[0]++;
			counts[1] ++;
		}
	}
}
