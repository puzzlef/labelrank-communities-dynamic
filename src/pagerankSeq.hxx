#pragma once
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "pagerank.hxx"

using std::vector;
using std::swap;




// PAGERANK-FACTOR
// ---------------
// For contribution factors of vertices (unchanging).

template <class T>
void pagerankFactorW(vector<T>& a, const vector<T>& vfact, const vector<int>& vdata, int N, T p) {
  for (int u=0; u<N; u++) {
    int d = vdata[u];
    T   f = vfact[u];
    a[u] = d>0? p*f/d : -p*f/N;
  }
}




// PAGERANK-TELEPORT
// -----------------
// For teleport contribution from vertices (inc. dead ends).

template <class T>
T pagerankTeleport(const vector<T>& r, const vector<T>& f, int N, T p) {
  T a = (1-p)/N;
  for (int u=0; u<N; u++)
    if (f[u]<0) a -= f[u]*r[u];
  return a;
}




// PAGERANK-CALCULATE
// ------------------
// For rank calculation from in-edges.

template <class T>
void pagerankCalculateW(vector<T>& a, const vector<T>& c, const vector<int>& vfrom, const vector<int>& efrom, int i, int n, T c0) {
  for (int v=i; v<i+n; v++)
    a[v] = c0 + sumValuesAt(c, sliceIterable(efrom, vfrom[v], vfrom[v+1]));
}




// PAGERANK-ERROR
// --------------
// For convergence check.

template <class T>
T pagerankError(const vector<T>& x, const vector<T>& y, int i, int N, int EF) {
  switch (EF) {
    case 1:  return l1Norm(x, y, i, N);
    case 2:  return l2Norm(x, y, i, N);
    default: return liNorm(x, y, i, N);
  }
}




// PAGERANK
// --------
// For Monolithic / Componentwise PageRank.

template <class H, class J, class M, class FL, class T=float>
PagerankResult<T> pagerankSeq(const H& xt, const J& ks, int i, const M& ns, FL fl, const vector<T> *q, const PagerankOptions<T>& o) {
  int  N  = xt.order();
  T    p  = o.damping;
  T    E  = o.tolerance;
  int  L  = o.maxIterations, l = 0;
  int  EF = o.toleranceNorm;
  auto vfrom = sourceOffsetsAs(xt, ks, int());
  auto efrom = destinationIndicesAs(xt, ks, int());
  auto vdata = vertexData(xt, ks);
  vector<T> a(N), r(N), c(N), f(N), qc, vfact;
  if (q) qc = compressContainer(xt, *q, ks);
  if (o.factors) vfact = compressContainer(xt, *o.factors, ks);
  else { vfact = vector<T>(N); fillValueU(vfact, T(1)); }
  float t = measureDurationMarked([&](auto mark) {
    if (q) copyValuesW(r, qc);    // copy old ranks (qc), if given
    else fillValueU(r, T(1)/N);
    copyValuesW(a, r);
    mark([&] { pagerankFactorW(f, vfact, vdata, N, p); multiplyValuesPositiveW(c, a, f, 0, N); });  // calculate factors (f) and contributions (c)
    mark([&] { l = fl(a, r, c, f, vfrom, efrom, i, ns, N, p, E, L, EF); });                         // calculate ranks of vertices
  }, o.repeat);
  return {decompressContainer(xt, a, ks), l, t};
}
