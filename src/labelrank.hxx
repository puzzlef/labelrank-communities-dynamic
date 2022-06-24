#pragma once




// LABELRANK-OPTIONS
// -----------------

template <class V>
struct LabelrankOptions {
  int repeat;
  int maxIterations;
  V   inflation;
  V   cutoff;
  V   conditionalUpdate;

  LabelrankOptions(int repeat, int maxIterations, V inflation, V cutoff, V conditionalUpdate) :
  repeat(repeat), maxIterations(maxIterations), inflation(inflation), cutoff(cutoff), conditionalUpdate(conditionalUpdate) {}
};




// LABELRANK-RESULT
// -----------------

template <class K>
struct LabelrankResult {
  vector<K> membership;
  int   iterations;
  float time;

  LabelrankResult(vector<K>&& membership, int iterations=0, float time=0) :
  membership(membership), iterations(iterations), time(time) {}

  LabelrankResult(vector<K>& membership, int iterations=0, float time=0) :
  membership(move(membership)), iterations(iterations), time(time) {}
};