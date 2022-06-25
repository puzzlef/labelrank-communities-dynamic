#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




template <class G>
void runExperiment(const G& x, int batch) {
  using T = float; vector<T> *init = nullptr;
  enum NormFunction { L0=0, L1=1, L2=2, Li=3 };
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int batch = argc>2? stoi(argv[2]) : 10;
  OutDiGraph<int, None, float> x;
  printf("Loading graph %s ...\n", file);
  readMtxW(x, file); println(x);
  auto fl = [](auto u) { return true; };
  selfLoopW(x, fl); print(x); printf(" (selfLoopAllVertices)\n");
  runExperiment(x, batch);
  printf("\n");
  return 0;
}
