

#include "Algebra.hpp"
#include "BasicDoFs.hpp"
#include "MatrixGraph.hpp"
#include "Preconditioner.hpp"
#include "Small.hpp"
#include "Solver.hpp"
#include "TestEnvironment.hpp"

std::pair<Matrix, Vector> assemble(int N) {
  int NN = (N - 2) * (N - 2);
  Meshes *m = new Meshes{"ST_squ_simple", 0, 0, 0};
  CellDoF *d = new CellDoF{NN};
  mout << (*m)[0].CellCount() << endl;
  MatrixGraph *mg = new MatrixGraph{(*m)[0], dof{d}};
  mout << (*mg).Size() << endl;
  VectorMatrixBase *mgg = new VectorMatrixBase{*mg};
  Matrix M{*mgg};
  double *a = M();
  mout << std::defaultfloat;
  Vector rhs{0.0, *mgg};
  double *rhs_ptr = rhs();
  for (int y = 0; y < NN; y++) {
    int vals = 0;
    for (int x = 0; x < NN; x++) {
      *a = 0;
      if (x == y) { *a = 4; }
      if ((x == y + 1 && x % (N - 2) > 0) || (x == y - 1 && x % (N - 2) != (N - 3))
          || x == y + N - 2 || x == y - N + 2) {
        *a = -1;
        vals++;
      }
      /*auto s= std::to_string(int(*a));
      s = s.size()==1?" "+s:s;
      mout << s << "  ";*/
      a++;
    }
    // mout << endl;
    *rhs_ptr = 4 - vals;
    rhs_ptr++;
  }
  return {M, rhs};
}

std::tuple<Matrix, Matrix, Matrix> partition(Matrix &M) {
  const int NN = M.size();
  Matrix L{M}, D{M}, U{M};
  L = 0;
  D = 0;
  U = 0;
  double *l = L(), *d = D(), *u = U(), *m = M();
  for (int y = 0; y < NN; y++) {
    for (int x = 0; x < NN; x++) {
      if (x < y) {
        *l = *m;
      } else if (x > y) {
        *u = *m;
      } else {
        *d = *m;
      }
      l++;
      u++;
      d++;
      m++;
    }
  }
  return {L, D, U};
}

void printMatrix(Matrix &M) {
  double *a = M();
  int NN = M.size();
  for (int y = 0; y < NN; y++) {
    for (int x = 0; x < NN; x++) {
      if (*a >= 0) {
        mout << " " << *a << " ";
      } else {
        mout << *a << " ";
      }
      a++;
    }
    mout << endl;
  }
}

class TestPreconditioner : public Test {

  void SetUp() override {
    int N = 10;
    auto [M, RHS] = assemble(N);
    int NN = RHS.size();

    auto PC = GetPC("Jacobi");
    PC->Construct(M);
    Vector x{4.0, RHS};
    Vector PCx = (*PC) * x;
    Vector PCx_exp{1.0, RHS};
    Vector difference{PCx - PCx_exp};
    ASSERT_LE(difference.norm(), 1e-14);

    auto [L, D, U] = partition(M);

    Matrix DL = D + L;
    Solver s;
    s(DL);
  }

  void TearDown() override {}
};

TEST_F(TestPreconditioner, TestPreconditioner) {}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}