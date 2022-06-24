#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;



// Find the best possible edge to insert (that reduces gini coefficient).
template <class G, class T=float>
auto edgeInsertBest(const G& x, const vector<T> *init=nullptr) {
  using K = typename G::key_type;
  T min = T(); K minu = K(), minv = K();
  auto a  = pagerankMonolithicSeq(x, init);
  auto lc = lorenzCurve(a.ranks);
  auto gc = giniCoefficient(lc);
  x.forEachVertexKey([&](auto u) {
    x.forEachVertexKey([&](auto v) {
      if (x.hasEdge(u, v)) return;
      auto y  = duplicate(x);
      y.addEdge(u, v);
      y.correct();
      auto b  = pagerankMonolithicSeq(y, &a.ranks);
      auto ld = lorenzCurve(b.ranks);
      auto gd = giniCoefficient(ld);
      if (gd>=gc || (minu && gd>=min)) return;
      min = gd; minu = u; minv = v;
    });
  });
  return make_pair(minu, minv);
}


// Process a batch with given edge adding process.
template <class G, class T, class F>
void runBatch(const G& x, const vector<T> *init, const char *name, int batch, F fn) {
  auto a  = pagerankMonolithicSeq(x, init);
  auto lc = lorenzCurve(a.ranks);
  auto gc = giniCoefficient(lc);
  printf("[%05d edges; %.10f gini_coef.] %s\n", 0, gc, name);
  auto y  = duplicate(x);
  auto r  = a.ranks;
  for (int i=0; i<batch; i++) {
    auto [u, v] = fn(y, r);
    y.addEdge(u, v);
    y.correct();
    auto b  = pagerankMonolithicSeq(y, &r);
    auto ld = lorenzCurve(b.ranks);
    auto gd = giniCoefficient(ld);
    printf("[%05d edges; %.10f gini_coef.] %s\n", i+1, gd, name);
    r = move(b.ranks);
  }
}


template <class G>
void runExperiment(const G& x, int batch) {
  using T = float; vector<T> *init = nullptr;
  enum NormFunction { L0=0, L1=1, L2=2, Li=3 };

  // Try different heuristics.
  const char *name;
  // name = "edgeInsertBest";
  // runBatch(x, init, name, batch, [&](auto y, auto r) {
  //   return edgeInsertBest(y, init);
  // });
  name = "edgeInsertCxrx";
  runBatch(x, init, name, batch, [&](auto y, auto r) {
    return edgeInsertCxrx(y, r);
  });
  name = "edgeInsertCxSx";
  runBatch(x, init, name, batch, [&](auto y, auto r) {
    auto yt = transpose(y);
    auto b  = pagerankMonolithicSeq(yt, init);
    return edgeInsertCxSx(y, r, b.ranks);
  });
  name = "edgeInsertCxSr";
  runBatch(x, init, name, batch, [&](auto y, auto r) {
    vector<T> f(y.span());
    y.forEachVertexKey([&](auto u) { f[u] = 1 - r[u]; });
    auto yt = transpose(y);
    auto b  = pagerankMonolithicSeq(yt, init, {&f});
    return edgeInsertCxSr(y, r, b.ranks);
  });
  name = "edgeInsertCRrx";
  runBatch(x, init, name, batch, [&](auto y, auto r) {
    return edgeInsertCRrx(y, r);
  });
  name = "edgeInsertCRSx";
  runBatch(x, init, name, batch, [&](auto y, auto r) {
    auto yt = transpose(y);
    auto b  = pagerankMonolithicSeq(yt, init);
    return edgeInsertCRSx(y, r, b.ranks);
  });
  name = "edgeInsertCRSr";
  runBatch(x, init, name, batch, [&](auto y, auto r) {
    vector<T> f(y.span());
    y.forEachVertexKey([&](auto u) { f[u] = 1 - r[u]; });
    auto yt = transpose(y);
    auto b  = pagerankMonolithicSeq(yt, init, {&f});
    return edgeInsertCRSr(y, r, b.ranks);
  });
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int batch = argc>2? stoi(argv[2]) : 10;
  printf("Loading graph %s ...\n", file);
  auto x  = readMtxOutDiGraph(file); println(x);
  selfLoopW(x, [](auto u) { return true; }); print(x); printf(" (selfLoopAllVertices)\n");
  runExperiment(x, batch);
  printf("\n");
  return 0;
}
