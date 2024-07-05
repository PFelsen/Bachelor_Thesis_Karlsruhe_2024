#include "KLExpansion.hpp"
#include <cmath>
#include "Assertion.hpp"

double KLExpansion::EvalSample(const Point &x) {
  double kl_evaluation = mean;
  for (size_t k = 0; k < sqrtEigenvalues.size(); k++) {
    kl_evaluation += sqrtEigenvalues[k] * Eigenfunction(x, k) * expansion_parameters[k];
  }
  return FieldTransform(kl_evaluation);
}

SimpleKLExpansion::SimpleKLExpansion() : KLExpansion() {
  mean = 0.0;
  sqrtEigenvalues = {1.0};
}

string SimpleKLExpansion::Name() const { return "Simple2DKLExpansion"; }

double SimpleKLExpansion::Eigenfunction(const Point &x, int k) { return 1; }

double SimpleKLExpansion::FieldTransform(double kl_evaluation) { return kl_evaluation; }

ExpKLExpansion::ExpKLExpansion() : KLExpansion() {
  mean = 0.0;
  sqrtEigenvalues = {1.0};
}

string ExpKLExpansion::Name() const { return "Exp2DKLExpansion"; }

double ExpKLExpansion::Eigenfunction(const Point &x, int k) { return 1; }

double ExpKLExpansion::FieldTransform(double kl_evaluation) { return std::exp(kl_evaluation); }

LinearKLExpansion::LinearKLExpansion() : KLExpansion() {
  mean = 1.1;
  sqrtEigenvalues = {1.0};
}

string LinearKLExpansion::Name() const { return "LinearExpansion"; }

double LinearKLExpansion::Eigenfunction(const Point &x, int k) { return 1; }

double LinearKLExpansion::FieldTransform(double kl_evaluation) { return kl_evaluation; }

TwoDimKLExpansion::TwoDimKLExpansion() : KLExpansion() {
  mean = 0.15;
  sqrtEigenvalues = {std::sqrt(1.2), std::sqrt(1.0), std::sqrt(0.9), std::sqrt(0.85)};
  expansion_parameters.resize(4);
}

string TwoDimKLExpansion::Name() const { return "TwoDimKLExpansion"; }

double TwoDimKLExpansion::Eigenfunction(const Point &x, int k) {
  double x0 = x[0];
  double x1 = x[1];
  switch (k) {
  case 0:
    return std::sin(M_PI * x0) * std::cos(M_PI * x1);
  case 1:
    return std::cos(M_PI * x0) * std::cos(M_PI * x1);
  case 2:
    return std::cos(M_PI * x0) * std::sin(M_PI * x1);
  case 3:
    return std::sin(M_PI * x0) * std::sin(M_PI * x1);
  default:
    return 0;
  }
}

double TwoDimKLExpansion::FieldTransform(double kl_evaluation) { return std::exp(kl_evaluation); }

ComplexTwoDimKLExpansion::ComplexTwoDimKLExpansion() : KLExpansion() {
  mean = 0.15;
  for (int j = 1; j <= 64; ++j) {
    sqrtEigenvalues.push_back(std::sqrt((2.0 / 5.0) * std::pow(4.0, -j)));
  }
  expansion_parameters.resize(64);
}

string ComplexTwoDimKLExpansion::Name() const { return "ComplexTwoDimKLExpansion"; }

double ComplexTwoDimKLExpansion::Eigenfunction(const Point &x, int k) {
  double x0 = x[0];
  double x1 = x[1];
  int i = k / 8;
  int i_f = i / 2 + 1;
  int j = k % 8;
  int j_f = j / 2 + 1;
  double f1 = (i % 2 == 1) ? std::sin(i_f * M_PI * x0) : std::cos(i_f * M_PI * x0);
  double f2 = (j % 2 == 1) ? std::sin(j_f * M_PI * x1) : std::cos(j_f * M_PI * x1);
  return f1 * f2;
}

double ComplexTwoDimKLExpansion::FieldTransform(double kl_evaluation) {
  return std::exp(kl_evaluation);
}

KLExpansion *CreateRandomField(const std::string &field_generator_name) {
  if (field_generator_name == "Simple2DKL") return new SimpleKLExpansion();
  if (field_generator_name == "Exp2DKL") return new ExpKLExpansion();
  if (field_generator_name == "Linear2DKL") return new LinearKLExpansion();
  if (field_generator_name == "Short2DKL") return new TwoDimKLExpansion();
  if (field_generator_name == "Long2DKL") return new ComplexTwoDimKLExpansion();
  Exit(field_generator_name + " not found");
}

std::unique_ptr<KLExpansion> CreateRandomFieldUnique(const std::string &field_generator_name) {
  return std::unique_ptr<KLExpansion>(CreateRandomField(field_generator_name));
}