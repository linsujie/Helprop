#ifndef NDINDEX_H
#define NDINDEX_H

#include<vector>

struct NDIndex {
  int dim;
  int n;
  std::vector<int> anti;
  std::vector<std::vector<std::vector<int> > > slice;
  std::vector<std::vector<int> > slice_anti;

  NDIndex() {};
  NDIndex(int dim);

  static inline bool is_up(int index, int ix) {
    return (index & (1 << ix));
  }

  void int2vec(int index, std::vector<bool>& vec) const;
  void show(int index) const;
};

#endif /* NDINDEX_H */