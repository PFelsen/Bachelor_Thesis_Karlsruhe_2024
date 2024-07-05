#ifndef RULE_H
#define RULE_H

#include "Celltype.hpp"

class Rule {
  const CELLTYPE tp;
public:
  const std::vector<short> node;
  std::vector<short> face;

  Rule(CELLTYPE TP, std::vector<short> cellNodes) : tp(TP), node(std::move(cellNodes)) {}

  Rule(const Rule &R) : tp(R.tp), node(R.node), face(R.face) {}

  Rule() : tp(NONE) {}

  CELLTYPE type() const { return tp; }

  template<typename T>
  void operator()(const std::vector<T> &z, std::vector<T> &x) const {
    x.resize(node.size());
    for (int j = 0; j < node.size(); ++j)
      x[j] = z[node[j]];
  }

  template<typename T>
  std::vector<Point> getCorners(const std::vector<T> &z) const {
    std::vector<Point> corners;
    (*this)(z, corners);
    return corners;
  }

  template<typename OSTREAM>
  OSTREAM &print(OSTREAM &os) const {
    os << "z ";
    for (int i = 0; i < node.size(); ++i)
      os << node[i] << " ";
    os << "f ";
    for (int i = 0; i < face.size(); ++i)
      os << face[i] << " ";
    return os << "tp " << tp;
  }

  int operator[](int i) const { return node[i]; }
};

inline std::ostream &operator<<(std::ostream &os, const Rule &R) { return R.print(os); }

#endif // RULE_H