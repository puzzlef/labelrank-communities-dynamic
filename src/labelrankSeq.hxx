#pragma once
#include <utility>
#include <cmath>
#include <vector>
#include "_main.hxx"
#include "labelrank.hxx"

using std::vector;
using std::pow;




/**
 * Get best label for a labelset,
 * @param x labelset
 * @return label with highest probability
 */
template <class B>
auto labelsetBestLabel(const B& x) {
  using K = typename B::key_type;
  using V = typename B::value_type;
  V max = V(); K maxk = K();
  x.forEach([&](auto k, auto v) {
    if (v<max) return;
    max = v; maxk = k;
  });
  return make_pair(maxk, max);
}


/**
 * Combine a labelset probabilities to target labelset with given weight.
 * @param a target labelset
 * @param x labelset to combine
 * @param w combining weight
 */
template <class B, class V>
void labelsetCombineU(B& a, const B& x, V w) {
  x.forEach([&](auto k, auto v) {
    if (a.has(k)) a.set(k, a.get(k) + w*v);
    else a.add(k, w*v);
  });
}


/**
 * Multiply a value with each probability in target labelset (probability scaling).
 * @param a target labelset
 * @param m value to multiply
 */
template <class B, class V>
void labelsetMultiplyU(B& a, V m) {
  a.forEachValue([&](auto& v) { v *= m; });
}


/**
 * Raise each probability with a given exponent in target labelset (inflation operator).
 * @param a target labelset
 * @param e exponent value
 */
template <class B, class V>
void labelsetPowU(B& a, V e) {
  a.forEachValue([&](auto& v) { v = pow(v, e); });
}


/**
 * Filter labels probabilities below given threshold in target labelset (cutoff operator).
 * @param a target labelset
 * @param th threshold value
 */
template <class B, class V>
void labelsetFilterAboveU(B& a, V th) {
  auto max = labelsetBestLabel(a).second, thmax = th*max;
  a.filterIfValue([&](auto v) { return v>=thmax; });
}


/**
 * Finish combining of labelsets to target labelset (probability scaling + inflation operator + cutoff operator).
 * @param a target labelset
 * @param m value to multiply
 * @param e exponent value
 * @param th threshold value
 */
template <class B, class V>
void labelsetCombineEndU(B& a, V m, V e, V th) {
  a.forEach([&](auto k, auto& v) { v = pow(v*m, e); });
  auto max = labelsetBestLabel(a).second, thmax = th*max;
  a.filterIfValue([&](auto v) { return v>=thmax; });
}


/**
 * Tells whether first labelset is subset of second.
 * @param x first labelset
 * @param y second labelset
 * @return is subset?
 */
template <class B>
bool labelsetIsSubset(const B& x, const B& y) {
  bool a = true;
  x.forEachKey([&](auto k) { a = a && y.has(k); });
  return a;
}




/**
 * Initialize labelset for a given vertex.
 * @param a target labelset
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 * @param th threshold value
 */
template <class B, class G, class K, class V>
void labelrankInitializeVertexW(B& a, const G& x, K u, V e, V th) {
  V sumw = V();
  x.forEachEdge(u, [&](auto v, auto w) {
    a.add(v, w);
    sumw += w;
  });
  labelsetCombineEndU(a, 1/sumw, e, th);
}


/**
 * Update labelset for a given vertex.
 * @param a target labelset
 * @param ls original labelsets
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 * @param th threshold value
 */
template <class B, class G, class K, class V>
void labelrankUpdateVertexW(B& a, const vector<B>& ls, const G& x, K u, V e, V th) {
  V sumw = V();
  x.forEachEdge(u, [&](auto v, auto w) {
    labelsetCombineU(a, ls[v], w);
    sumw += w;
  });
  labelsetCombineEndU(a, 1/sumw, e, th);
}


/**
 * Check if a vertex is stable.
 * @param ls labelsets
 * @param x original graph
 * @param u given vertex
 * @param q conditional update parameter
 * @return true
 * @return false
 */
template <class B, class G, class K, class V>
bool labelrankIsVertexStable(const B& ls, const G& x, K u, V q) {
  K count = K();
  x.forEachEdgeKey(u, [&](auto v) {
    if (labelsetIsSubset(ls[u], ls[v])) count++;
  });
  return count > q * x.degree(u);
}




template <class G, class V=float>
auto labelrankSeq(const G& x, const LabelrankOptions<V>& o={}) {
  using K = typename G::key_type;
  vector<OrderedBitset<K, V>> ls(x.span());
  vector<OrderedBitset<K, V>> ms(x.span());
  x.forEachVertexKey([&](auto u) {
    labelrankInitializeVertexW(ls, x, u, o.inflation, o.cutoff);
  });
  int i = 0;
  K updatedPrev = K();
  while (true) {
    K updated = K();
    size_t labels = 0;
    labelrankClearVertices(ms, x);
    x.forEachVertexKey([&](auto u) {
      labels += ls[u].size();
      // if (labelrankIsVertexStable(ls, x, u, o.conditionalUpdate)) ms[u] = ls[u];
      // else { labelrankUpdateVertexW(ms, ls, x, u, o.inflation, o.cutoff); updated++; }
      labelrankUpdateVertexW(ms, ls, x, u, o.inflation, o.cutoff); updated++;
      if (u%1000==0) printf("update vertex: %zu, labels: %zu\n", u, labels);
    }); i++;
    swap(ls, ms);
    printf("i: %d, updated: %d\n", i, updated);
    if (!updated || updated==updatedPrev) break;
    updatedPrev = updated;
  }
  vector<K> a(x.span()); K zeros = K();
  x.forEachVertexKey([&](auto u) {
    a[u] = labelsetBestLabel(ls[u]).first;
    if (a[u]==0) zeros++;
  });
  printf("zeros: %d\n", zeros);
  return LabelrankResult(a, i, 0);
}
