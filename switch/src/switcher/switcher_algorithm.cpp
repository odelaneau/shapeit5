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

#include <switcher/switcher_header.h>

#include <models/mendel_solver.h>
#include <models/genotype_checker.h>
#include <models/haplotype_checker.h>

using namespace std;

void switcher::process() {
	mendel_solver MP(H);
	genotype_checker GC(H);
	haplotype_checker HC(H, options["nbins"].as < int > ());

	//
	if (options.count("pedigree")) {
		MP.solve(options.count("singleton"));
		MP.count();
		MP.writePerSample(options["output"].as < string > () + ".sample.mendel.txt.gz");
		MP.writePerVariant(options["output"].as < string > () + ".variant.mendel.txt.gz");
		MP.writeImbalance(options["output"].as < string > () + ".variant.imbalance.txt.gz");
		MP.writePedigree(options["output"].as < string > () + ".sample.pedigree");
		//MP.writeDistances(options["output"].as < string > () + ".sample.distance.txt.gz");
	} else MP.set();

	//
	GC.check();
	GC.writePerSample(options["output"].as < string > () + ".sample.typing.txt.gz");
	GC.writePerVariant(options["output"].as < string > () + ".variant.typing.txt.gz");

	//
	HC.check();
	HC.writePerSample(options["output"].as < string > () + ".sample.switch.txt.gz");
	HC.writePerVariant(options["output"].as < string > () + ".variant.switch.txt.gz");
	HC.writePerFrequency(options["output"].as < string > () + ".frequency.switch.txt.gz");
	HC.writePerType(options["output"].as < string > () + ".type.switch.txt.gz");
	HC.writeFlipSwitchErrorPerSample(options["output"].as < string > () + ".flipsAndSwitches.txt.gz");
	HC.writeBlock(options["output"].as < string > () + ".block.switch.txt.gz");
	HC.writeCalibration(options["output"].as < string > () + ".calibration.switch.txt.gz");

}
