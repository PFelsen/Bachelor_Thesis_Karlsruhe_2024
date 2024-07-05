#include "DebuggingTools.hpp"

void printFirstBlock(const Matrix &M) {

  const Scalar *m = M.GetData()();
  mout << std::fixed;

  for (int i = 0; i < 1; ++i) {
    int n = M.Dof(i);
    for (int k = 0; k < n; ++k) {
      for (int l = 0; l < n; ++l, ++m) {
        mout << *m << " ";
      }
      mout << endl;
    }
  }
  mout << std::defaultfloat;
}