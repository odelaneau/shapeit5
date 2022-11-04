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

#include <io/sparse2plain.h>
#include <io/plain2sparse.h>

using namespace std;

void converter::convert() {

	int mode1 = options.count("input-plain") + options.count("output-sparse");
	int mode2 = options.count("input-sparse") + options.count("output-plain");

	if (mode1==2 && mode2==0) {
		plain2sparse(options["input-plain"].as < string > (), options["output-sparse"].as < string > (), options["region"].as < string > (), options["thread"].as < int > (), options["maf"].as < double > ()).convert();
	} else if (mode1==0 && mode2==2) {
		sparse2plain(options["output-plain"].as < string > (), options["input-sparse"].as < string > (), options["region"].as < string > (), options["thread"].as < int > ()).convert();
	}

}
