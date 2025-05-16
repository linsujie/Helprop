#include<vector>
#include<iostream>
#include "NDIndex.h"

using namespace std;

NDIndex::NDIndex(int dim) : dim(dim), n(1 << dim) {
  anti.resize(n);
  for (int index = 0; index < n; index++)
    anti[index] = n - 1 - index;
  
  slice.resize(dim);
  for (int ix = 0; ix < dim; ix++) {
    slice[ix].resize(2);
    for (int i = 0; i < n; i++)
      is_up(i, ix) ? slice[ix][1].push_back(i) : slice[ix][0].push_back(i);
  }

  slice_anti.resize(dim);
  for (int ix = 0; ix < dim; ix++) {
    slice_anti[ix].resize(n);
    for (int i = 0; i < n; i++)
      if (is_up(i, ix))
        slice_anti[ix][i] = i - (1 << ix);
      else
        slice_anti[ix][i] = i + (1 << ix);
  }
}

void NDIndex::int2vec(int index, std::vector<bool>& vec) const {
  vec.resize(dim);
  for (int ix = 0; ix < vec.size(); ix++)
    vec[ix] = is_up(index, ix);
}

void NDIndex::show(int index) const {
  vector<bool> vec(dim);
  int2vec(index, vec);

  for (int ix = 0; ix < vec.size(); ix++)
    cout << vec[ix];
}
