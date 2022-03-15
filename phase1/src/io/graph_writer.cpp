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
#include <io/graph_writer.h>

graph_writer::graph_writer(genotype_set & _G, variant_map & _V): G(_G), V(_V) {
}

graph_writer::~graph_writer() {
}

void graph_writer::binary_write(output_file & fout, const vector<bool> & x) {
	vector<bool>::size_type n = x.size();
	fout.write((const char*)&n, sizeof(std::vector<bool>::size_type));
	for(std::vector<bool>::size_type i = 0; i < n;) {
		unsigned char aggr = 0;
		for(unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1) if(x.at(i)) aggr |= mask;
		fout.write((const char*)&aggr, sizeof(unsigned char));
    }
}

void graph_writer::string_write(output_file & fout, string & x) {
	size_t size_str = x.size();
	fout.write(reinterpret_cast<char*>(&size_str), sizeof(size_str));
	fout.write(reinterpret_cast<char*>(&x[0]), size_str);
}

void graph_writer::writeGraphs(string fname) {
	// Init
	tac.clock();
	output_file fd (fname);

	//Write variant map
	int n_variants = V.vec_pos.size();
	fd.write(reinterpret_cast<char*>(&n_variants), sizeof(n_variants));
	for (int l = 0 ; l < n_variants ; l ++) {
		string_write(fd, V.vec_pos[l]->chr);
		fd.write(reinterpret_cast<char*>(&V.vec_pos[l]->bp), sizeof(V.vec_pos[l]->bp));
		string_write(fd, V.vec_pos[l]->id);
		string_write(fd, V.vec_pos[l]->ref);
		string_write(fd, V.vec_pos[l]->alt);
		fd.write(reinterpret_cast<char*>(&V.vec_pos[l]->idx), sizeof(V.vec_pos[l]->idx));
	}

	//Write genotype graphs
	fd.write(reinterpret_cast<char*>(&G.n_ind), sizeof(G.n_ind));
	for (int g  = 0 ; g < G.n_ind ; g++) {
		// name
		string_write(fd, G.vecG[g]->name);
		// integers
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->index), sizeof(G.vecG[g]->index));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->n_segments), sizeof(G.vecG[g]->n_segments));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->n_variants), sizeof(G.vecG[g]->n_variants));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->n_ambiguous), sizeof(G.vecG[g]->n_ambiguous));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->n_missing), sizeof(G.vecG[g]->n_missing));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->n_transitions), sizeof(G.vecG[g]->n_transitions));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->n_stored_transitionProbs), sizeof(G.vecG[g]->n_stored_transitionProbs));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->n_storage_events), sizeof(G.vecG[g]->n_storage_events));

		// vectors
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->Variants[0]), G.vecG[g]->Variants.size());
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->Ambiguous[0]), G.vecG[g]->Ambiguous.size());
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->Diplotypes[0]), G.vecG[g]->Diplotypes.size() * sizeof(unsigned long));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->Lengths[0]), G.vecG[g]->Lengths.size() * sizeof(unsigned short));
		binary_write(fd, G.vecG[g]->ProbMask);
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->ProbStored[0]), G.vecG[g]->ProbStored.size() * sizeof(float));
		fd.write(reinterpret_cast<char*>(&G.vecG[g]->ProbMissing[0]), G.vecG[g]->ProbMissing.size() * sizeof(float));
	}
	vrb.bullet("BIN writing [Compressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)");
}
