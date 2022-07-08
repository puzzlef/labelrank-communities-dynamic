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
template <class BA, class BX, class V>
void labelsetCombineU(BA& a, const BX& x, V w) {
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
 * Finish combining of labelsets to target labelset (probability scaling + inflation operator + cutoff operator).
 * @param a target labelset
 * @param m value to multiply
 * @param e exponent value
 * @param th threshold value
 */
template <class B, class V>
void labelsetCombineEndU2(B& a, V m, V e, V th) {
  // auto t0 = timeNow();
  a.forEach([&](auto k, auto& v) { v = v*m; });  // v = pow(v*m, e);
  // auto t1 = timeNow();
  auto max = labelsetBestLabel(a).second, thmax = th*max;
  // auto t2 = timeNow();
  a.filterIfValue([&](auto v) { return v>=thmax; });
  // auto t3 = timeNow();
  // auto d0 = durationMilliseconds(t0, t1);
  // auto d1 = durationMilliseconds(t1, t2);
  // auto d2 = durationMilliseconds(t2, t3);
  // printf("-> d0: %fms, d1: %fms, d2: %fms\n", d0, d1, d2);
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
 * @param a temporary labelset
 * @param as target labelsets
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 * @param th threshold value
 */
template <class BT, class B, class G, class K, class V>
void labelrankInitializeVertexW(BT& a, vector<B>& as, const G& x, K u, V e, V th) {
  V sumw = V(); a.clear();
  x.forEachEdge(u, [&](auto v, auto w) {
    a.add(v, w);
    sumw += w;
  });
  labelsetCombineEndU(a, 1/sumw, e, th);
  copyW(as[u], a);
}


/**
 * Update labelset for a given vertex.
 * @param a temporary labelset
 * @param as target labelsets
 * @param ls original labelsets
 * @param x original graph
 * @param u given vertex
 * @param e exponent value
 * @param th threshold value
 */
template <class BT, class B, class G, class K, class V>
void labelrankUpdateVertexW(BT& a, vector<B>& as, const vector<B>& ls, const G& x, K u, V e, V th) {
  auto t0 = timeNow();
  V sumw = V(); a.clear();
  x.forEachEdge(u, [&](auto v, auto w) {
    labelsetCombineU(a, ls[v], w);
    sumw += w;
  });
  auto t1 = timeNow();
  labelsetCombineEndU2(a, 1/sumw, e, th);
  auto t2 = timeNow();
  copyW(as[u], a);
  auto t3 = timeNow();
  auto d0 = durationMilliseconds(t0, t1);
  auto d1 = durationMilliseconds(t1, t2);
  auto d2 = durationMilliseconds(t2, t3);
  printf("d0: %fms, d1: %fms, d2: %fms\n", d0, d1, d2);
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
  DenseBitset<K, V> la(x.span());
  vector<UnorderedBitset<K, V>> ls(x.span());
  vector<UnorderedBitset<K, V>> ms(x.span());
  auto start = timeNow();
  x.forEachVertexKey([&](auto u) {
    labelrankInitializeVertexW(la, ls, x, u, o.inflation, o.cutoff);
  });
  auto stop = timeNow();
  printf("init_time: %fms\n", durationMilliseconds(start, stop));
  int i = 0;
  K updatedPrev = K();
  while (true) {
    K updated = K();
    auto start = timeNow();
    x.forEachVertexKey([&](auto u) {
      if (labelrankIsVertexStable(ls, x, u, o.conditionalUpdate)) ms[u] = ls[u];
      else { labelrankUpdateVertexW(la, ms, ls, x, u, o.inflation, o.cutoff); updated++; }
    }); i++;
    auto stop = timeNow();
    swap(ls, ms);
    printf("i: %d, updated: %d, time: %fms\n", i, updated, durationMilliseconds(start, stop));
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
