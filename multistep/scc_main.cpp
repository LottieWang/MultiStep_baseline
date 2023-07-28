/*
//@HEADER
// *****************************************************************************
//
//       Multistep: (Strongly) Connected Components Algorithms
//              Copyright (2016) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  George M. Slota (gmslota@sandia.gov)
//                      Siva Rajamanickam (srajama@sandia.gov)
//                      Kamesh Madduri (madduri@cse.psu.edu)
//
// *****************************************************************************
//@HEADER
*/


using namespace std;

#include <cstdlib>
#include <fstream>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <omp.h>
#include <fstream>
#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/random.h>
// #include "pbbslib/sample_sort.h"

#define VERBOSE 0
#define DEBUG 0
#define VERIFY 0
#define TIMING 1
#define TRIM_LEVEL 1

#define THREAD_QUEUE_SIZE 2048
#define ALPHA 15.0
#define BETA 25


double timer()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}

typedef struct graph {
  size_t n, m;
  int* out_array;
  int* in_array;
  unsigned* out_degree_list;
  unsigned* in_degree_list;
} graph;
#define out_degree(g, n) (g.out_degree_list[n+1] - g.out_degree_list[n])
#define in_degree(g, n) (g.in_degree_list[n+1] - g.in_degree_list[n])
#define out_vertices(g, n) &g.out_array[g.out_degree_list[n]]
#define in_vertices(g, n) &g.in_array[g.in_degree_list[n]]

#include "scc_trim.cpp"
#include "scc_fwbw.cpp"
#include "scc_color.cpp"
#include "scc_serial.cpp"
#include "scc_verify.cpp"

void read_edge(char* filename,
  int& n, unsigned& m,
  int*& srcs, int*& dsts)
{
  ifstream infile;
  string line;
  infile.open(filename);

  getline(infile, line, ' ');
  n = atoi(line.c_str());
  getline(infile, line);
  m = strtoul(line.c_str(), NULL, 10);

  int src, dst;
  unsigned counter = 0;

  srcs = new int[m];
  dsts = new int[m];
  for (unsigned i = 0; i < m; ++i)
  {
    getline(infile, line, ' ');
    src = atoi(line.c_str());
    getline(infile, line);
    dst = atoi(line.c_str());

    srcs[counter] = src;
    dsts[counter] = dst;
    ++counter;
  }

  infile.close();
}

void read_binary(char const *filename, size_t &n, size_t &m, int* &out_array, int* &in_array, unsigned* &out_degree_list, unsigned* &in_degree_list, int &max_deg_vert, double &avg_degree) {
  std::ifstream ifs(filename);
  if(!ifs.is_open()) {
    std::cerr << "Error: file " << filename << " does not exist\n";
    exit(EXIT_FAILURE);
  }
  size_t _n, _m, sizes;
  ifs.read(reinterpret_cast<char*>(&_n), sizeof(size_t));
  ifs.read(reinterpret_cast<char*>(&_m), sizeof(size_t));
  ifs.read(reinterpret_cast<char*>(&sizes), sizeof(size_t));
  n = _n, m = _m;
  assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);

  parlay::sequence<uint64_t> offset(n + 1);
  parlay::sequence<uint32_t> es(m);
  ifs.read(reinterpret_cast<char*>(offset.begin()), (n + 1) * 8);
  ifs.read(reinterpret_cast<char*>(es.begin()), m * 4);
  if(ifs.peek() != EOF) {
    fprintf(stderr, "Error: Bad data\n");
    exit(EXIT_FAILURE);
  }
  ifs.close();

  // out
  out_degree_list = new unsigned[n+1];
  out_array = new int[m];
  parlay::parallel_for(0, n+1, [&](size_t i) {
    out_degree_list[i] = offset[i];
  });
  parlay::parallel_for(0, m, [&](size_t i) {
    out_array[i] = es[i];
  });

  // in
  in_degree_list = new unsigned[n+1];
  in_array = new int[m];

  //for(int i = 0; i <= n; i++) {
    //in_degree_list[i] = 0;
  //}
  //for(int i = 0; i < n; i++) {
    //for(size_t j = offset[i]; j < offset[i+1]; j++) {
      //in_degree_list[es[j]]++;
    //}
  //}
  //uint32_t sum = 0;
  //for(int i = 0; i < n; i++) {
    //uint32_t tmp = in_degree_list[i];
    //in_degree_list[i] = sum;
    //sum += tmp;
  //}
  //in_degree_list[n] = sum;
  //assert(sum == m);
  //unsigned* tmp = new unsigned[n];
  //copy(in_degree_list, in_degree_list+n, tmp);
  //for(int i = 0; i < n; i++) {
    //for(size_t j = offset[i]; j < offset[i+1]; j++) {
      //in_array[tmp[es[j]]++] = i;
    //}
  //}

  parlay::sequence<std::pair<unsigned, unsigned>> edge_list(m);
  parlay::parallel_for(0, n, [&](size_t i) {
    parlay::parallel_for(offset[i], offset[i+1], [&](size_t j) {
      edge_list[j] = {es[j], i};
    });
  });
  // sample_sort_inplace(edge_list.slice(), [&](std::pair<unsigned, unsigned> a, std::pair<unsigned, unsigned> b) {
  //   return a < b;
  // });
  parlay::integer_sort_inplace(edge_list,
																	[&](const std::pair<unsigned, unsigned>& p) { return p.first; });
  parlay::parallel_for(0, edge_list[0].first, [&](size_t i) {
    in_degree_list[i] = 0;
  });
  in_array[0]=edge_list[0].second;
  parlay::parallel_for(1, m, [&](size_t i) {
    unsigned u = edge_list[i].first;
    unsigned v = edge_list[i].second;
    in_array[i] = v;
    if(edge_list[i-1].first != u) {
      in_degree_list[u] = i;
    }
    if(i == m-1 || edge_list[i+1].first != u) {
      int end = (i==m-1?n+1:edge_list[i+1].first);
      parlay::parallel_for(u+1, end, [&](size_t j) {
        in_degree_list[j] = i+1;
      });
    }
  });
  
  max_deg_vert = 0;
  avg_degree = 0;
  double max_deg = 0;
  for(size_t i = 0; i < n; i++) {
    int in_deg = in_degree_list[i+1]-in_degree_list[i];
    int out_deg = out_degree_list[i+1]-out_degree_list[i];
    double prod = (double)in_deg*out_deg;
    if(prod > max_deg) {
      max_deg = prod;
      max_deg_vert = i;
    }
    avg_degree += out_deg;
  }
  avg_degree /= n;


  // validate
  //parallel_for(0, n+1, [&](size_t i) {
    //assert(i == n || in_degree_list[i] <= in_degree_list[i+1]);
    //assert(in_degree_list[i] <= m);
    //assert(i == n || out_degree_list[i] <= out_degree_list[i+1]);
    //assert(out_degree_list[i] <= m);
  //});
  //parallel_for(0, m, [&](size_t i) {
    //assert(in_array[i] < n);
    //assert(out_array[i] < n);
  //});
  //parallel_for(0, n, [&](size_t i) {
    //parallel_for(offset[i], offset[i+1], [&](size_t j) {
      //assert(j == offset[i+1]-1 || es[j] <= es[j+1]);
    //});
  //});
};

void create_csr(int n, unsigned m, 
  int* srcs, int* dsts,
  int*& out_array, int*& in_array,
  unsigned*& out_degree_list, unsigned*& in_degree_list,
  int& max_deg_vert, double& avg_degree)
{
  out_array = new int[m];
  in_array = new int[m];
  out_degree_list = new unsigned[n+1];
  in_degree_list = new unsigned[n+1];

  for (unsigned i = 0; i < m; ++i)
    out_array[i] = 0;
  for (unsigned i = 0; i < m; ++i)
    in_array[i] = 0;  
  for (int i = 0; i < n+1; ++i)
    out_degree_list[i] = 0;
  for (int i = 0; i < n+1; ++i)
    in_degree_list[i] = 0;

  unsigned* temp_counts = new unsigned[n];
  for (int i = 0; i < n; ++i)
    temp_counts[i] = 0;
  for (unsigned i = 0; i < m; ++i)
    ++temp_counts[srcs[i]];
  for (int i = 0; i < n; ++i)
    out_degree_list[i+1] = out_degree_list[i] + temp_counts[i];
  copy(out_degree_list, out_degree_list + n, temp_counts);
  for (unsigned i = 0; i < m; ++i)
    out_array[temp_counts[srcs[i]]++] = dsts[i];

  for (int i = 0; i < n; ++i)
    temp_counts[i] = 0;
  for (unsigned i = 0; i < m; ++i)
    ++temp_counts[dsts[i]];
  for (int i = 0; i < n; ++i)
    in_degree_list[i+1] = in_degree_list[i] + temp_counts[i];
  copy(in_degree_list, in_degree_list + n, temp_counts);
  for (unsigned i = 0; i < m; ++i)
    in_array[temp_counts[dsts[i]]++] = srcs[i];
  delete [] temp_counts;

  avg_degree = 0.0;
  double max_degree = 0.0;
  for (int i = 0; i < n; ++i)
  {
    unsigned out_degree = out_degree_list[i+1] - out_degree_list[i];
    unsigned in_degree = in_degree_list[i+1] - in_degree_list[i];
    double degree = (double)out_degree * (double)in_degree;
    avg_degree += (double)out_degree;
    if (degree > max_degree)
    {
      max_deg_vert = i;
      max_degree = degree;
    }
  }
  avg_degree /= (double)n;

#if DEBUG
  int max_out = out_degree_list[max_deg_vert+1] - out_degree_list[max_deg_vert];
  int max_in = in_degree_list[max_deg_vert+1] - in_degree_list[max_deg_vert];
  int max_tot = max_out > max_in ? max_out : max_in;
  printf("max deg vert: %d (%d), avg_degree %9.2lf\n", max_deg_vert, max_tot, avg_degree);
#endif
}

void output_scc(graph& g, int* scc_maps, char* output_file)
{
  std::ofstream outfile;
  outfile.open(output_file);

  for (size_t i = 0; i < g.n; ++i)
    outfile << scc_maps[i] << std::endl;

  outfile.close();
}

void print_usage(char** argv)
{
  printf("Usage: %s [graph] [optional: output file]\n", argv[0]);
  exit(0);
}

void run_scc(graph& g, int*& scc_maps,
  int max_deg_vert, double avg_degree, int vert_cutoff)
{

  int num_trim = 0;
  int num_fwbw = 0;
  int num_color = 0;
  int num_serial = 0;

#if VERBOSE
  printf("Performing Trim step ... \n");
  double elt = timer();
  double start_time = elt;
#endif

  scc_maps = new int[g.n];
  bool* valid = new bool[g.n];
  int* valid_verts = new int[g.n];
  int num_valid = g.n;

#pragma omp parallel
{
#pragma omp for nowait
  for (size_t i = 0; i < g.n; ++i) valid[i] = true;
#pragma omp for nowait
  for (size_t i = 0; i < g.n; ++i) scc_maps[i] = -1;
#pragma omp for nowait
  for (size_t i = 0; i < g.n; ++i) valid_verts[i] = i;
}


if (TRIM_LEVEL == 0)
{
/*  num_trim = scc_trim_none(g, valid,
    valid_verts, num_valid,
    scc_maps);*/
}
else if (TRIM_LEVEL == 1)
{
  num_trim = scc_trim(g, valid, 
    valid_verts, num_valid, 
    scc_maps);
}
else if (TRIM_LEVEL == 2)
{
  num_trim = scc_trim_complete(g, valid,
    valid_verts, num_valid,
    scc_maps);
  num_valid = 0;
  for (size_t i = 0; i < g.n; ++i) if (valid[i]) valid_verts[num_valid++] = i;
}

#if VERIFY
  scc_verify(g, scc_maps);  
#endif

#if VERBOSE
  elt = timer() - elt;
  printf("\tDone, %9.6lf, %d verts trimmed\n", elt, num_trim);
  printf("Performing FWBW step ... \n");
  elt = timer();
#endif

  num_fwbw = scc_fwbw(g, valid, 
    valid_verts, num_valid, 
    max_deg_vert, avg_degree,
    scc_maps);

#if VERIFY
  scc_verify(g, scc_maps);  
#endif

#if VERBOSE
  elt = timer() - elt;
  printf("\tDone, %9.6lf, %d verts found\n", elt, num_fwbw);
  printf("Performing Coloring step ... \n");
  elt = timer();
#endif
  
  num_color = scc_color(g, valid,
    valid_verts, num_valid, 
    vert_cutoff,
    scc_maps);

#if VERIFY
  scc_verify(g, scc_maps);  
#endif

#if VERBOSE
  elt = timer() - elt;
  printf("\tDone, %9.6lf, %d sccs found\n", elt, num_color);
  printf("Performing Serial step ... \n");
  elt = timer();
#endif
 
  tarjan t(g, valid,
    valid_verts, num_valid, 
    scc_maps);
  num_serial = t.run_valid();

#if VERBOSE
  elt = timer() - elt;
  printf("\tDone, %9.6lf, %d sccs found\n", elt, num_serial);
  start_time = timer() - start_time;
  printf("TOTAL: %9.6lf\n", start_time);
#endif

}

int main(int argc, char** argv)
{
  setbuf(stdout, NULL);
  if (argc < 2)
    print_usage(argv);
  int rounds = 10;
  bool stats = false;
  char* out_file = NULL;
  for (int i = 2; i< argc; i++){
    string s(argv[i]);
    if (s=="-t"){
      rounds = atoi(argv[++i]);
    }else if (s=="-stats"){
      stats=true;
    }else if (s=="-o"){
      out_file = argv[++i];
    }
  }
  //int* srcs;
  //int* dsts;
  size_t n, m;
  int* out_array;
  int* in_array;
  unsigned* out_degree_list;
  unsigned* in_degree_list;
  int* scc_maps;
  int max_deg_vert;
  double avg_degree;
  int vert_cutoff = 10000;

#if VERBOSE
  double elt, start_time;
  printf("Reading %s ... ", argv[1]);
  elt = timer();
#endif

  //read_edge(argv[1], n, m, srcs, dsts);

#if VERBOSE
  elt = timer() - elt;
  printf("Done, %9.6lf\n", elt);
  printf("Creating graph ... ");
  elt = timer();
#endif

  read_binary(argv[1], n, m, out_array, in_array, out_degree_list, in_degree_list, max_deg_vert, avg_degree);
  std::cout << "Graph: " << argv[1] << " n: " << n << " m: " << m << std::endl;
  //create_csr(n, m, srcs, dsts, 
    //out_array, in_array,
    //out_degree_list, in_degree_list,
    //max_deg_vert, avg_degree);
  graph g = {n, m, out_array, in_array, out_degree_list, in_degree_list};
  //delete [] srcs;
  //delete [] dsts;

#if VERBOSE
  elt = timer() - elt;
  printf("Done, %9.6lf, n: %d, m: %u\n", elt, n, m);
  printf("Doing Multistep ... ");
  elt = timer();
#endif

#if TIMING
  double exec_time = timer();
#endif

  double temp_time = timer();
  run_scc(g, scc_maps,
    max_deg_vert, avg_degree, vert_cutoff);
  temp_time = timer() - temp_time;
  printf("round 0: %f\n", temp_time);
  double total_time = 0;
  for(int i = 0; i < rounds; i++) {
    double exec_time = timer();
    run_scc(g, scc_maps, max_deg_vert, avg_degree, vert_cutoff);
    exec_time = timer() - exec_time;
    printf("round %d: %f\n", i+1, exec_time);
    total_time += exec_time;
  }
  printf("average cost: %f\n", total_time/rounds);

#if TIMING
  exec_time = timer() - exec_time;
  //printf("Multistep SCC time: %9.6lf\n", exec_time);
#endif
if (stats){
  scc_verify(g, scc_maps);
}

#if VERBOSE
  elt = timer() - elt;
  printf("Done, %9.6lf\n", elt);
#endif

  if (out_file!= NULL)
    output_scc(g, scc_maps, out_file);

  delete [] out_array;
  delete [] out_degree_list;
  delete [] scc_maps;

  return 0;
}
