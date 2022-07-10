#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




template <class G>
void runExperiment(const G& x, int repeat) {
  using K = typename G::key_type;
  auto M = edgeWeight(x)/2;
  printf("[original_modularity: %f]\n", modularity(x, M, 1.0f));
  LabelrankResult<K> a = labelrankSeq<4>(x);
  printf("[%09.3f ms; %03d iters.] labelrankSeq\n", a.time, a.iterations);
  auto fc = [&](auto u) { return a.membership[u]; };
  printf("[modularity: %f]\n", modularity(x, fc, M, 1.0f));
}


int main(int argc, char **argv) {
  using K = int;
  using V = float;
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  OutDiGraph<K, None, V> x; V w = 1;
  printf("Loading graph %s ...\n", file);
  readMtxW(x, file); println(x);
  auto y  = symmetricize(x); print(y); printf(" (symmetricize)\n");
  auto fl = [](auto u) { return true; };
  selfLoopU(y, w, fl); print(y); printf(" (selfLoopAllVertices)\n");
  runExperiment(y, repeat);
  printf("\n");
  return 0;
}
