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

#include <converter/converter_header.h>

#include <modes/bcf2binary.h>
#include <modes/binary2bcf.h>

using namespace std;

void converter::convert() {
	//Retrieve parameter values
	string region = options["region"].as < string > ();
	string format = options["format"].as < string > ();
	string finput = options["input"].as < string > ();
	string foutput = options["output"].as < string > ();
	uint32_t nthreads = options["thread"].as < int > ();
	float maf = options["maf"].as < float > ();

	//
	if (format == "bg") bcf2binary (region, maf, nthreads, CONV_BCF_BG).convert(finput, foutput);

	//
	else if (format == "bh") bcf2binary (region, maf, nthreads, CONV_BCF_BH).convert(finput, foutput);

	//
	else if (format == "sg") bcf2binary (region, maf, nthreads, CONV_BCF_SG).convert(finput, foutput);

	//
	else if (format == "sh") bcf2binary (region, maf, nthreads, CONV_BCF_SH).convert(finput, foutput);

	//
	else  if (isBCF(format)) binary2bcf (region, nthreads).convert(finput, foutput);

	//
	else vrb.error("Output format [" + format + "] unrecognized");
}
