#pragma once
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
  return maxk;
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
  using K = typename B::key_type; vector<K> removes;
  a.forEach([&](auto k, auto v) { if (v<th) removes.push_back(k); });
  for (auto k : removes)
    a.remove(k);
  // TODO: bitset.filter()
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
  using K = typename B::key_type; vector<K> removes;
  a.forEach([&](auto k, auto& v) {
    v = pow(v*m, e);
    if (v<th) removes.push_back(k);
  });
  for (auto k : removes)
    a.remove(k);
  // TODO: bitset.filter()
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
 * @param as target labelsets
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 * @param th threshold value
 */
template <class B, class G, class K, class V>
void labelrankInitializeVertexW(vector<B>& as, const G& x, K u, V e, V th) {
  V sumw = V();
  x.forEachEdge(u, [&](auto v, auto w) {
    as[u].add(v, w);
    sumw += w;
  });
  labelsetCombineEndU(as[u], 1/sumw, e, th);
}


/**
 * Update labelset for a given vertex.
 * @param as target labelsets
 * @param ls original labelsets
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 * @param th threshold value
 */
template <class B, class G, class K, class V>
void labelrankUpdateVertexW(vector<B>& as, const vector<B>& ls, const G& x, K u, V e, V th) {
  V sumw = V();
  x.forEachEdge(u, [&](auto v, auto w) {
    labelsetCombineU(as[u], ls[v], w);
    sumw += w;
  });
  labelsetCombineEndU(as[u], 1/sumw, e, th);
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


/**
 * Clear labelsets for all vertices.
 * @param as target labelsets
 * @param x original graph
 */
template <class B, class G>
void labelrankClearVertices(B& as, const G& x) {
  x.forEachVertexKey([&](auto u) { as[u].clear(); });
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
  K oldupdated = K();
  while (true) {
    K updated = K();
    labelrankClearVertices(ms, x);
    x.forEachVertexKey([&](auto u) {
      if (labelrankIsVertexStable(ls, x, u, o.conditionalUpdate)) ms[u] = ls[u];
      else { labelrankUpdateVertexW(ms, ls, x, u, o.inflation, o.cutoff); updated++; }
    }); i++;
    swap(ls, ms);
    printf("i: %d, updated: %d\n", i, updated);
    if (!updated || updated==oldupdated) break;
    oldupdated = updated;
  }
  vector<K> a(x.span()); K zeros = K();
  x.forEachVertexKey([&](auto u) {
    a[u] = labelsetBestLabel(ls[u]);
    if (a[u]==0) zeros++;
  });
  printf("zeros: %d\n", zeros);
  return LabelrankResult(a, i, 0);
}
