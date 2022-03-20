////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <models/hmm_scaffold/hmm_scaffold_header.h>

/*
 * ALL THESE ROUTINES CAN BE VECTORIZED BIG TIME ...
 */

void hmm_scaffold::forward_initTransitions(unsigned int h) {
	fill(alpha.begin() + sstates[h], alpha.begin() + sstates[h] + nstates[h], 1.0f / nstates[h]);
}


void hmm_scaffold::forward_updateTransitions(unsigned int vs, unsigned int h) {
	assert(vs > 0);
	float f0 = M.t[vs-1] / nstates[h];
	float f1 = M.nt[vs-1];
	for (int k = 0 ; k < nstates[h] ; k ++) alpha[sstates[h]+k] = alpha[sstates[h]+k] * f0 + f1;
}

void hmm_scaffold::forward_updateEmission(unsigned int vs, unsigned int h) {
	for (int k = 0 ; k < nstates[h] ; k ++)
		alpha[sstates[h]+k] *= emit[bufferA0[h]!=bufferA0[K[sstates[h]+k]]];
}

void hmm_scaffold::forward_normalize(unsigned int vs, unsigned int h) {
	float scale = 0.0f;
	for (int k = 0 ; k < nstates[h] ; k ++) scale += alpha[sstates[h]+k];
	scale = 1.0f / scale;
	for (int k = 0 ; k < nstates[h] ; k ++) alpha[sstates[h]+k] *= scale;
}

void hmm_scaffold::forward_reverseEmission(unsigned int vs, unsigned int h) {
	for (int k = 0 ; k < nstates[h] ; k ++)
		alpha[sstates[h]+k] *= rev_emit[bufferA0[h]!=bufferA0[K[sstates[h]+k]]];
}

void hmm_scaffold::forward_reverseTransitions(unsigned int vs, unsigned int h) {
	float f0 = nstates[h] / M.t[vs];
	float f1 = M.nt[vs];
	for (int k = 0 ; k < nstates[h] ; k ++)
		alpha[sstates[h]+k] = (alpha[sstates[h]+k]  - f1) * f0;
}

void hmm_scaffold::backward_initTransitions(unsigned int h) {
	fill(beta.begin() + sstates[h], beta.begin() + sstates[h] + nstates[h], 1.0f / nstates[h]);
}

void hmm_scaffold::backward_updateTransitions(unsigned int vs, unsigned int h) {
	assert(vs < (C.n_scaffold_variants-1));
	float f0 = M.t[vs] / nstates[h];
	float f1 = M.nt[vs];
	for (int k = 0 ; k < nstates[h] ; k ++) beta[sstates[h]+k] = beta[sstates[h]+k] * f0 + f1;
}

void hmm_scaffold::backward_updateEmission(unsigned int vs, unsigned int h) {
	for (int k = 0 ; k < nstates[h] ; k ++)
		beta[sstates[h]+k] *= emit[bufferA0[h]!=bufferA0[K[sstates[h]+k]]];
}

void hmm_scaffold::backward_normalize(unsigned int vs, unsigned int h) {
	float scale = 0.0f;
	for (int k = 0 ; k < nstates[h] ; k ++) scale += beta[sstates[h]+k];
	scale = 1.0f / scale;
	for (int k = 0 ; k < nstates[h] ; k ++) beta[sstates[h]+k] *= scale;
}

void hmm_scaffold::getAlphaBetaProduct(unsigned int h, vector < float > & alphaXbeta) {
	float scale = 0.0f;
	for (int k = 0 ; k < nstates[h] ; k ++) {
		//cout << h << " " << k << " " << alpha[sstates[h]+k] << " " << beta[sstates[h]+k] << endl;
		alphaXbeta[k] = alpha[sstates[h]+k] * beta[sstates[h]+k];
		scale += alphaXbeta[k];
	}
	scale = 1.0f / scale;
	for (int k = 0 ; k < nstates[h] ; k ++) alphaXbeta[k] *= scale;
}



