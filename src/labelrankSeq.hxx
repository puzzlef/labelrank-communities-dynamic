#pragma once
#include <utility>
#include <cmath>
#include <array>
#include "_main.hxx"
#include "Labelset.hxx"
#include "labelrank.hxx"

using std::pair;
using std::array;
using std::pow;




/**
 * Initialize labelset for a given vertex.
 * @param a accumulator labelset (scratch)
 * @param as target labelsets
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 */
template <class G, class K, class V, size_t N>
void labelrankInitializeVertexW(ALabelset<K, V>& a, vector<Labelset<K, V, N>>& as, const G& x, K u, V e) {
  V sumw = V(); a.clear();
  x.forEachEdge(u, [&](auto v, auto w) {
    a.set(v, w);
    sumw += w;
  });
  labelsetReorderU(a);
  labelsetCopyW(as[u], a);
  labelsetMultiplyPowU(as[u], 1/sumw, e);
}


/**
 * Update labelset for a given vertex.
 * @param a accumulator labelset (scratch)
 * @param as target labelsets
 * @param ls original labelsets
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 */
template <class G, class K, class V, size_t N>
void labelrankUpdateVertexW(ALabelset<K, V>& a, vector<Labelset<K, V, N>>& as, const vector<Labelset<K, V, N>>& ls, const G& x, K u, V e) {
  V sumw = V(); a.clear();
  x.forEachEdge(u, [&](auto v, auto w) {
    labelsetCombineU(a, ls[v], w);
    sumw += w;
  });
  labelsetReorderU(a);
  labelsetCopyW(as[u], a);
  labelsetMultiplyPowU(as[u], 1/sumw, e);
}


/**
 * Check if a vertex is stable.
 * @param ls labelsets
 * @param x original graph
 * @param u given vertex
 * @param q conditional update parameter
 */
template <class B, class G, class K, class V>
bool labelrankIsVertexStable(const B& ls, const G& x, K u, V q) {
  K count = K();
  x.forEachEdgeKey(u, [&](auto v) {
    if (labelsetIsSubset(ls[u], ls[v])) count++;
  });
  return count > q * x.degree(u);
}




template <size_t N, class G>
auto labelrankSeq(const G& x, const LabelrankOptions& o={}) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  ALabelset<K, V> la(x.span());
  vector<Labelset<K, V, N>> ls(x.span());
  vector<Labelset<K, V, N>> ms(x.span());
  auto t0 = timeNow();
  x.forEachVertexKey([&](auto u) {
    labelrankInitializeVertexW(la, ls, x, u, V(o.inflation));
  });
  auto t1 = timeNow();
  auto d0 = durationMilliseconds(t0, t1);
  printf("init_time: %fms\n", d0);
  int i = 0;
  size_t updatedPrev = 0;
  while (true) {
    size_t updated = 0;
    auto t2 = timeNow();
    x.forEachVertexKey([&](auto u) {
      if (labelrankIsVertexStable(ls, x, u, o.conditionalUpdate)) ms[u] = ls[u];
      else { labelrankUpdateVertexW(la, ms, ls, x, u, o.inflation); updated++; }
    }); i++;
    auto t3 = timeNow();
    auto d2 = durationMilliseconds(t2, t3);
    printf("i: %d, updated: %zu, time: %fms\n", i, updated, d2);
    swap(ls, ms);
    if (!updated || updated==updatedPrev) break;
    updatedPrev = updated;
  }
  vector<K> a(x.span()); size_t zeros = 0;
  x.forEachVertexKey([&](auto u) {
    a[u] = ls[u][0].first;
    if (a[u]==0) ++zeros;
  });
  printf("zeros: %zu\n", zeros);
  return LabelrankResult(a, i, 0);
}
