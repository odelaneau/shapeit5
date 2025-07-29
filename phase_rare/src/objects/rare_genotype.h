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

#ifndef _RARE_GENOTYPE_H
#define _RARE_GENOTYPE_H

#include <utils/otools.h>
#include <utils/sparse_genotype.h>

#define SETBIT(n,i)	(n)|=(1U<<(i));
#define CLRBIT(n,i)	(n)&=~(1U<<i);
#define GETBIT(n,i)	(((n)>>(i))&1U);


class rare_genotype : public sparse_genotype {
public:

	static float ee;
	static float ed;

	rare_genotype(uint32_t _idx, bool _het, bool _mis, bool _al0, bool _al1, bool _pha) : sparse_genotype(_idx, _het,  _mis, _al0, _al1, _pha) {
		if (!pha && al0 != al1) {
			if (rng.flipCoin()) { al0 = 0; al1 = 1; }
			else { al0 = 1; al1 = 0; }
		}
	}

	rare_genotype(uint32_t _val) : sparse_genotype(_val) {
	}	
	
	void randomize() {
		if (!pha) {
			if (het) sparse_genotype::phase(rng.getInt(2) + 1);
			if (mis) sparse_genotype::phase(rng.getInt(4));
		}
	}

	void phase(uint32_t g) {
		sparse_genotype::phase(g);
	}

	void phase(float prb0, float prb1) {
		if (!pha && (het || mis)) {
			std::vector < double > gprobs = std::vector < double > (4, 0.0f);
			float p01 = std::max(prb0, std::numeric_limits<float>::min());
			float p00 = std::max(1.0f - prb0, std::numeric_limits<float>::min());
			float p11 = std::max(prb1, std::numeric_limits<float>::min());
			float p10 = std::max(1.0f - prb1, std::numeric_limits<float>::min());

			if (het) {
				gprobs[0] = 0.0f;
				gprobs[1] = (p00*ee + p01*ed) * (p10*ed + p11*ee);
				gprobs[2] = (p00*ed + p01*ee) * (p10*ee + p11*ed);
				gprobs[3] = 0.0f;
			} else {
				gprobs[0] = (p00*ee + p01*ed) * (p10*ee + p11*ed);
				gprobs[1] = (p00*ee + p01*ed) * (p10*ed + p11*ee);
				gprobs[2] = (p00*ed + p01*ee) * (p10*ee + p11*ed);
				gprobs[3] = (p00*ed + p01*ee) * (p10*ed + p11*ee);
			}

			int32_t maxg = alg.imax(gprobs);
			switch (maxg) {
			case 0:	al0 = 0; al1 = 0; break;
			case 1:	al0 = 0; al1 = 1; break;
			case 2:	al0 = 1; al1 = 0; break;
			case 3:	al0 = 1; al1 = 1; break;
			}
			prob = gprobs[maxg] / (gprobs[0] + gprobs[1] + gprobs[2] + gprobs[3]);
		}
	}

	void impute(float prb0, float prb1) {
		if (mis) {
			std::vector < double > gprobs = std::vector < double > (4, 0.0f);
			float p01 = std::max(prb0, std::numeric_limits<float>::min());
			float p00 = std::max(1.0f - prb0, std::numeric_limits<float>::min());
			float p11 = std::max(prb1, std::numeric_limits<float>::min());
			float p10 = std::max(1.0f - prb1, std::numeric_limits<float>::min());

			gprobs[0] = (p00*ee + p01*ed) * (p10*ee + p11*ed);
			gprobs[1] = 0.0f;
			gprobs[2] = 0.0f;
			gprobs[3] = (p00*ed + p01*ee) * (p10*ed + p11*ee);

			int32_t maxg = alg.imax(gprobs);
			switch (maxg) {
			case 0:	al0 = 0; al1 = 0; break;
			case 1:	al0 = 0; al1 = 1; break;
			case 2:	al0 = 1; al1 = 0; break;
			case 3:	al0 = 1; al1 = 1; break;
			}
			prob = gprobs[maxg] / (gprobs[0] + gprobs[1] + gprobs[2] + gprobs[3]);
		}
	}

	std::string str() {
		return std::string("[" + stb.str(idx) + "/" + stb.str(het) + "/" + stb.str(mis) + "/" + stb.str(al0) + "/" + stb.str(al1) + "/" + stb.str(pha) + "]");
	}
};

#endif
