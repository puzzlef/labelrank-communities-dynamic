#pragma once
#include <vector>
#include <utility>

using std::vector;
using std::make_pair;




// Pick vertex with low vertex property.
template <class G, class F>
auto vertexLowBy(const G& x, F fn) {
  using K = typename G::key_type;
  using T = decltype(fn(K()));
  T min  = T();
  K minu = K();
  x.forEachVertexKey([&](auto u) {
    T val = fn(u);
    if (val<0) return;
    if (minu && val>=min) return;
    min = val; minu = u;
  });
  return minu;
}


// Pick vertex with high vertex property.
template <class G, class F>
auto vertexHighBy(const G& x, F fn) {
  using K = typename G::key_type;
  using T = decltype(fn(K()));
  T max  = T();
  K maxu = K();
  x.forEachVertexKey([&](auto u) {
    T val = fn(u);
    if (val<0) return;
    if (maxu && val<=max) return;
    max = val; maxu = u;
  });
  return maxu;
}


template <class G, class K, class F>
auto inedgeHighBy(const G& x, K v, F fn) {
  using T = decltype(fn(K()));
  T max  = T();
  K maxu = K();
  x.forEachInEdgeKey(v, [&](auto u) {
    T val = fn(u);
    if (val<0) return;
    if (maxu && val<=max) return;
    max = val; maxu = u;
  });
  return maxu;
}




// Pick vertex with low rank.
template <class G, class T, class FT>
inline auto vertexLowRank(const G& x, const vector<T>& r, FT ft) {
  auto fn = [&](auto u) { return ft(u)? r[u] : -1; };
  return vertexLowBy(x, fn);
}
template <class G, class T>
inline auto vertexLowRank(const G& x, const vector<T>& r) {
  auto ft = [](auto u) { return true; };
  return vertexLowRank(x, r, ft);
}


// Pick vertex with high rank.
template <class G, class T, class FT>
inline auto vertexHighRank(const G& x, const vector<T>& r, FT ft) {
  auto fn = [&](auto u) { return ft(u)? r[u] : -1; };
  return vertexHighBy(x, fn);
}
template <class G, class T>
inline auto vertexHighRank(const G& x, const vector<T>& r) {
  auto ft = [](auto u) { return true; };
  return vertexHighRank(x, r, ft);
}


// Pick vertex with high contrib.
template <class G, class T, class FT>
inline auto vertexHighContrib(const G& x, const vector<T>& r, FT ft) {
  auto fn = [&](auto u) { return ft(u)? r[u] / (x.degree(u) + 1) : -1; };
  return vertexHighBy(x, fn);
}
template <class G, class T>
inline auto vertexHighContrib(const G& x, const vector<T>& r) {
  auto ft = [](auto u) { return true; };
  return vertexHighContrib(x, r, ft);
}


// Pick vertex with high contrib to high rank.
template <class G, class T, class FT>
inline auto vertexHighContribToHighRank(const G& x, const vector<T>& r, FT ft) {
  auto v  = vertexHighRank(x, r);
  auto fn = [&](auto u) { return ft(u)? r[u] / (x.degree(u) + 1) : -1; };
  return inedgeHighBy(x, v, fn);
}
template <class G, class T>
inline auto vertexHighContribToHighRank(const G& x, const vector<T>& r) {
  auto ft = [](auto u) { return true; };
  return vertexHighContribToHighRank(x, r, ft);
}


// Generate insertion edge from high contrib to low rank.
template <class G, class T>
auto edgeInsertCxrx(const G& x, const vector<T>& r) {
  auto u  = vertexHighContrib(x, r);
  auto ft = [&](auto v) { return !x.hasEdge(u, v); };
  auto v  = vertexLowRank(x, r, ft);
  return make_pair(u, v);
}

// Generate insertion edge from high contrib to high reverse rank.
template <class G, class T>
auto edgeInsertCxSx(const G& x, const vector<T>& r, const vector<T>& s) {
  auto u  = vertexHighContrib(x, r);
  auto ft = [&](auto v) { return !x.hasEdge(u, v); };
  auto v  = vertexHighRank(x, s, ft);
  return make_pair(u, v);
}

// Generate insertion edge from high contrib to high reverse rank with low ranks.
template <class G, class T>
inline auto edgeInsertCxSr(const G& x, const vector<T>& r, const vector<T>& s) {
  return edgeInsertCxSx(x, r, s);
}

// Generate insertion edge from high contrib targeting high rank to low rank.
template <class G, class T>
auto edgeInsertCRrx(const G& x, const vector<T>& r) {
  auto u  = vertexHighContribToHighRank(x, r);
  auto ft = [&](auto v) { return !x.hasEdge(u, v); };
  auto v  = vertexLowRank(x, r, ft);
  return make_pair(u, v);
}

// Generate insertion edge from high contrib targeting high rank to high reverse rank.
template <class G, class T>
auto edgeInsertCRSx(const G& x, const vector<T>& r, const vector<T>& s) {
  auto u  = vertexHighContribToHighRank(x, r);
  auto ft = [&](auto v) { return !x.hasEdge(u, v); };
  auto v  = vertexHighRank(x, s, ft);
  return make_pair(u, v);
}

// Generate insertion edge from high contrib targeting high rank to high reverse rank with low ranks.
template <class G, class T>
inline auto edgeInsertCRSr(const G& x, const vector<T>& r, const vector<T>& s) {
  return edgeInsertCRSx(x, r, s);
}
