#ifndef RANDOM_FIELD_HPP
#define RANDOM_FIELD_HPP

#include <string>
#include <vector>
#include "Point.hpp"
#include "RVector.hpp"
#include "SymmetricCovariance.hpp"

class KLExpansion {
protected:
  std::vector<double> sqrtEigenvalues;

  RVector expansion_parameters;

  std::string fieldType;

  double mean;
public:
  KLExpansion() : expansion_parameters(1) {}

  virtual double Eigenfunction(const Point &x, int k) = 0;

  virtual double FieldTransform(double kl_evaluation) = 0;

  virtual double EvalSample(const Point &x);

  void SetParameterVector(RVector parameters) { expansion_parameters = parameters; }

  virtual RVector &GetParameterVector() { return expansion_parameters; }

  virtual string Name() const = 0;
};

class SimpleKLExpansion : public KLExpansion {
public:
  SimpleKLExpansion();

  double Eigenfunction(const Point &x, int k) override;

  double FieldTransform(double kl_evaluation) override;

  string Name() const override;
};

class ExpKLExpansion : public KLExpansion {
public:
  ExpKLExpansion();

  double Eigenfunction(const Point &x, int k) override;

  double FieldTransform(double kl_evaluation) override;

  string Name() const override;
};

class LinearKLExpansion : public KLExpansion {
public:
  LinearKLExpansion();

  double Eigenfunction(const Point &x, int k) override;

  double FieldTransform(double kl_evaluation) override;

  string Name() const override;
};

class TwoDimKLExpansion : public KLExpansion {
public:
  TwoDimKLExpansion();

  double Eigenfunction(const Point &x, int k) override;

  double FieldTransform(double kl_evaluation) override;

  string Name() const override;
};

class ComplexTwoDimKLExpansion : public KLExpansion {
public:
  ComplexTwoDimKLExpansion();

  double Eigenfunction(const Point &x, int k) override;

  double FieldTransform(double kl_evaluation) override;

  string Name() const override;
};

KLExpansion *CreateRandomField(const std::string &field_generator_name);

std::unique_ptr<KLExpansion> CreateRandomFieldUnique(const std::string &field_generator_name);

#endif // RANDOM_FIELD_HPP
