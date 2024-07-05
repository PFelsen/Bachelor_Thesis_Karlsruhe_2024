#ifndef TESTTIMESTEPPING_HPP
#define TESTTIMESTEPPING_HPP

#include "Assemble.hpp"
#include "LagrangeDiscretization.hpp"
#include "MeshesCreator.hpp"
#include "TimeIntegrator.hpp"

class TestAssemble : public ILinearTimeAssemble {
public:
  TestAssemble() : ILinearTimeAssemble(0.0, 1.0, 1) {}

  vector<double> t_vec;
  vector<Vector> u_vec;

  void SetInitialValue(Vector &u) override {
    t_vec.push_back(Time());
    u_vec.push_back(u);
  }

  const char *Name() const override { return "TestAssemble"; };

  void Initialize(Vector &u) const override {}

  // Identity Matrix
  void MassMatrix(Matrix &massMatrix) const override {
    massMatrix = 0.0;
    for (row r = massMatrix.rows(); r != massMatrix.rows_end(); ++r)
      massMatrix(r, r)[0] = 1;
  }

  // Stiffness Matrix with -1 entries
  void SystemMatrix(Matrix &systemMatrix) const override {
    systemMatrix = 0.0;
    for (row r = systemMatrix.rows(); r != systemMatrix.rows_end(); ++r)
      systemMatrix(r, r)[0] = -1;
  }

  void RHS(double t, Vector &rhs) const override { rhs.Clear(); }

  void FinishTimeStep(const Vector &u) override {
    t_vec.push_back(Time());
    u_vec.push_back(u);
  }

  virtual void Solution(double t, Vector &uEx) {
    Scalar c = exp(-t);
    uEx *= c;
  };
};

template<int degree>
class TestPAssemble : public TestAssemble {
public:
  TestPAssemble() : TestAssemble() {}

  vector<double> t_vec;
  vector<Vector> u_vec;

  void SetInitialValue(Vector &u) override {
    t_vec.push_back(Time());
    u_vec.push_back(u);
  }

  // Identity Matrix
  void MassMatrix(Matrix &massMatrix) const override {
    massMatrix = 0.0;
    for (row r = massMatrix.rows(); r != massMatrix.rows_end(); ++r)
      massMatrix(r, r)[0] = 1;
  }

  void Solution(double t, Vector &uEx) {
    uEx = u_vec[0];
    Vector uTemp(1.0 / degree * pow(t, degree), uEx);
    uTemp.MakeAdditive();
    uEx += uTemp;
  };

  void RHS(double t, Vector &rhs) const override {
    rhs = degree > 0 ? pow(t, degree - 1) : 0.0;
    rhs.MakeAdditive();
  };

  void SystemMatrix(Matrix &systemMatrix) const override { systemMatrix = 0.0; }

  void FinishTimeStep(const Vector &u) override {
    t_vec.push_back(Time());
    u_vec.push_back(u);
  }
};

#endif // TESTTIMESTEPPING_HPP
