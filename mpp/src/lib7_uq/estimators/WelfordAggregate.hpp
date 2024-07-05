#ifndef WELFORDAGGREGATE_HPP
#define WELFORDAGGREGATE_HPP

#include <utility>

#include "Sample.hpp"
#include "Parallel.hpp"


struct SampleCounter {

  int M = 0;

  int Mpara = 0;

  int Mcomm = 0;

  friend std::ostream &operator<<(std::ostream &s, const SampleCounter &ctr) {
    return s << "M=" << ctr.M
             << " Mpara=" << ctr.Mpara
             << " Mcomm=" << ctr.Mcomm;
  }
};

class WelfordAggregate;

struct Errors {

  double mse = 0.0;

  double rmse = 0.0;

  double disc = 0.0;

  double disc2 = 0.0;

  double input = 0.0;

  void EstimateErrors(const WelfordAggregate &aggregate);

  static double EstimateNumeric(const WelfordAggregate &aggregate);

  static double EstimateStochastic(const WelfordAggregate &aggregate);

  friend std::ostream &operator<<(std::ostream &s, const Errors &errors) {
    return s << "rmse=" << errors.rmse
             << " mse=" << errors.mse
             << " input=" << errors.input
             << " disc^2=" << errors.disc2;
  }
};

template<typename T>
struct WelfordDelta {
  T one;
  T two;
  T three;
  T four;

  WelfordDelta(const T &newValue);

  WelfordDelta(const T &newValue, const T &oldMean);

  void ComputeDeltas();

  friend std::ostream &operator<<(std::ostream &s, const WelfordDelta &welfordDelta) {
    return s << "delta1=" << welfordDelta.one
             << " delta2=" << welfordDelta.two
             << " delta3=" << welfordDelta.three
             << " delta4=" << welfordDelta.four
             << endl;
  }
};

template<typename T>
struct WelfordData {
  T mean;
  T sVar;
  T skew;
  T kurt;
  T powSum2;
  T powSum3;
  T powSum4;
  double weight;
  double weight2;

  WelfordData(const T &newValue);

  WelfordData(const WelfordData<double> &welfordData);

  WelfordData(double value, const WelfordData<double> &welfordData);

  WelfordData(double value, const WelfordData<Vector> &welfordData);

  WelfordData(double value, const WelfordData<Vector> &welfordData,
              int commSplit);

  WelfordData(T mean, T sVar, T skew, T kurt, T powSum2, T powSum3, T powSum4,
              double weight, double weight2);

  T computeSqrt(const T &value) const;

  inline void ResetPowerSums() {
    mean = 0.0;
    powSum2 = 0.0;
    powSum3 = 0.0;
    powSum4 = 0.0;
  }

  inline void computeMean(double lW, double rW, const WelfordDelta<T> &delta) {
    double Wlr = (lW + rW);
    T deltaTmp = delta.one;
    deltaTmp *= (rW / Wlr);
    mean += deltaTmp;
  }

  inline void computeSVar() {
    sVar = powSum2;
    sVar *= 1.0 / (weight - weight / weight2);
  }

  inline void computeSkew() {
    skew = powSum3;
    skew *= sqrt(weight);
    T denominator = computeSqrt(powSum2);
    denominator *= powSum2;
    skew /= denominator;
  }

  inline void computeKurt() {
    kurt = powSum4;
    kurt *= weight;
    kurt /= powSum2;
    kurt /= powSum2;
  }

  inline void computePowSum2(double lW, double rW, const WelfordDelta<T> &delta) {
    double Wlr = (lW + rW);
    T delta2Tmp = delta.two;
    delta2Tmp *= (lW * rW) / Wlr;
    powSum2 += delta2Tmp;
  }

  inline void computePowSum2(double lW, double rW,
                             const WelfordData<T> &otherData,
                             const WelfordDelta<T> &delta) {
    computePowSum2(lW, rW, delta);
    powSum2 += otherData.powSum2;
  }

  inline void computePowSum3(double lW, double rW,
                             const WelfordDelta<T> &delta) {
    double Wlr = (rW + lW);
    T ldelta3Tmp = delta.three;
    T rdelta3Tmp = delta.three;
    ldelta3Tmp *= lW * rW * rW * (-rW) / (Wlr * Wlr * Wlr);
    rdelta3Tmp *= rW * lW * lW * lW / (Wlr * Wlr * Wlr);

    T ldeltaTmp = delta.one;
    ldeltaTmp *= 3 * (-rW) / Wlr;
    ldeltaTmp *= powSum2;

    powSum3 += ldelta3Tmp;
    powSum3 += rdelta3Tmp;
    powSum3 += ldeltaTmp;
  }

  inline void computePowSum3(double lW, double rW,
                             const WelfordData<T> &otherData,
                             const WelfordDelta<T> &delta) {
    computePowSum2(lW, rW, delta);
    double Wlr = (rW + lW);
    T rdeltaTmp = delta.one;
    rdeltaTmp *= 3 * lW / Wlr;
    rdeltaTmp *= otherData.powSum2;

    powSum3 += otherData.powSum3;
    powSum3 += rdeltaTmp;
  }

  inline void computePowSum4(double lW, double rW,
                             const WelfordDelta<T> &delta) {
    double Wlr = (rW + lW);
    T ldelta4Tmp = delta.four;
    T rdelta4Tmp = delta.four;
    ldelta4Tmp *= lW * pow(rW / Wlr, 4);
    rdelta4Tmp *= rW * pow(lW / Wlr, 4);

    T ldeltaTmp = delta.one;
    ldeltaTmp *= 4 * (-rW) / Wlr;
    ldeltaTmp *= powSum3;

    T ldelta2Tmp = delta.two;
    ldelta2Tmp *= 6 * pow(rW / Wlr, 2);
    ldelta2Tmp *= powSum2;

    powSum4 += ldelta4Tmp;
    powSum4 += rdelta4Tmp;
    powSum4 += ldelta2Tmp;
    powSum4 += ldeltaTmp;
  }

  inline void computePowSum4(double lW, double rW,
                             const WelfordData<T> &otherData,
                             const WelfordDelta<T> &delta) {
    computePowSum4(lW, rW, delta);
    double Wlr = (rW + lW);
    T rdeltaTmp = delta.one;
    rdeltaTmp *= 4 * lW / Wlr;
    rdeltaTmp *= otherData.powSum3;

    T rdelta2Tmp = delta.two;
    rdelta2Tmp *= 6 * pow(rW / Wlr, 2);
    rdelta2Tmp *= otherData.powSum2;

    powSum4 += otherData.powSum4;
    powSum4 += rdeltaTmp;
    powSum4 += rdelta2Tmp;
  }


  friend std::ostream &operator<<(std::ostream &s, const WelfordData &welfordData) {
    return s << "mean=" << welfordData.mean
             << " sVar=" << welfordData.sVar
             << " skew=" << welfordData.skew
             << " kurt=" << welfordData.kurt
             << endl;
  }
};

template<typename T>
class Aggregate {

  std::shared_ptr<WelfordData<T>> comm;

  std::shared_ptr<WelfordData<T>> para;

  std::shared_ptr<WelfordData<T>> total;

public:
  Aggregate() = default;

  Aggregate(const WelfordData<double> &total) :
      total(std::make_shared<WelfordData<double>>(total)) {}

  T SumOnCommSplit(const T &value, int commSplit);

  void UpdateTemplateOnComm(const T &newValue, double newWeight);

  void UpdateParallel(int maxCommSplit);

  void UpdateTotal();

  const WelfordData<T> &GetTotal() const {
    if (total) return *total;
    Exit("Total WelfordData not initialized.")
  }

  const WelfordData<T> &GetPara() const {
    if (para) return *para;
    Exit("Para WelfordData not initialized.")
  }

  const WelfordData<T> &GetComm() const {
    if (comm) return *comm;
    Exit("Comm WelfordData not initialized.")
  }

  T GetMean() const {
    if (total) return (total->mean);
    if constexpr (std::is_same<T, double>::value) return 0.0;
    else Exit("Total not initialized")
  }

  T GetSVar() const {
    if (total) return (total->sVar);
    if constexpr (std::is_same<T, double>::value) return 0.0;
    else Exit("Total not initialized")
  }

  T GetSkew() const {
    if (total) return (total->skew);
    if constexpr (std::is_same<T, double>::value) return 0.0;
    else Exit("Total not initialized")
  }

  T GetKurt() const {
    if (total) return (total->kurt);
    if constexpr (std::is_same<T, double>::value) return 0.0;
    else Exit("Total not initialized")
  }

  T GetPowSum2() const {
    if (total) return (total->powSum2);
    if constexpr (std::is_same<T, double>::value) return 0.0;
    else Exit("Total not initialized")
  }

  T GetPowSum3() const {
    if (total) return (total->powSum3);
    if constexpr (std::is_same<T, double>::value) return 0.0;
    else Exit("Total not initialized")
  }

  T GetPowSum4() const {
    if (total) return (total->powSum4);
    if constexpr (std::is_same<T, double>::value) return 0.0;
    else Exit("Total not initialized")
  }

  bool TotalInitialized() const {
    if (total) return true;
    return false;
  }

  bool ParaInitialized() const {
    if (para) return true;
    return false;
  }

  bool CommInitialized() const {
    if (comm) return true;
    return false;
  }

  double GetWeight() const { return total->weight; }

  double GetWeight2() const { return total->weight2; }

  friend WelfordAggregate;

};

struct QuantityOfInterest : public Aggregate<double> {

  QuantityOfInterest() = default;

  QuantityOfInterest(const WelfordData<double> &total) : Aggregate(total) {}

  void UpdateOnComm(const SampleSolution &fSol, const SampleSolution &cSol) {
    UpdateTemplateOnComm(fSol.Q, fSol.W);
  }
};

struct DeltaQuantityOfInterest : public Aggregate<double> {

  DeltaQuantityOfInterest() = default;

  DeltaQuantityOfInterest(const WelfordData<double> &total) : Aggregate(total) {}

  void UpdateOnComm(const SampleSolution &fSol, const SampleSolution &cSol) {
    UpdateTemplateOnComm((fSol.Q - cSol.Q), fSol.W);
  }
};

struct CostInSeconds : public Aggregate<double> {

  CostInSeconds() = default;

  CostInSeconds(const WelfordData<double> &total) : Aggregate(total) {}

  void UpdateOnComm(const SampleSolution &fSol, const SampleSolution &cSol) {
    UpdateTemplateOnComm((fSol.C + cSol.C), 1.0);
  }
};

#ifdef AGGREGATE_FOR_SOLUTION
#include "Transfers.hpp"

struct AggregateSolution : public Aggregate<Vector> {
  void UpdateOnComm(const SampleSolution &fineSample, const SampleSolution &coarseSample) {
    UpdateTemplateOnComm(fineSample.solution.vector, fineSample.W);
  }
};

struct DeltaSolution : public Aggregate<Vector> {
  void UpdateOnComm(const SampleSolution &fineSample, const SampleSolution &coarseSample) {
    if (fineSample.solution.vector.Level() == coarseSample.solution.vector.Level()) {
      UpdateTemplateOnComm(fineSample.solution.vector, fineSample.W);
      return;
    }
    Vector prolongatedCoarseSampleWithoutBC(fineSample.solution.vector);
    Vector fineSampleWithoutBC(fineSample.solution.vector);
    prolongatedCoarseSampleWithoutBC.ClearDirichletFlags();
    fineSampleWithoutBC.ClearDirichletFlags();
    auto transfer = GetTransfer(coarseSample.solution.vector, prolongatedCoarseSampleWithoutBC, "Matrix");
    transfer->Prolongate(coarseSample.solution.vector, prolongatedCoarseSampleWithoutBC);
    prolongatedCoarseSampleWithoutBC *= -1;
    prolongatedCoarseSampleWithoutBC += fineSampleWithoutBC;
    UpdateTemplateOnComm(prolongatedCoarseSampleWithoutBC, fineSample.W);
  }
};
#endif

class WelfordAggregate {
public:
  SampleCounter ctr;

  Errors errors;

  CostInSeconds C;

#ifdef AGGREGATE_FOR_SOLUTION
  DeltaSolution V;

  AggregateSolution U;
#endif

  QuantityOfInterest Q;

  DeltaQuantityOfInterest Y;

  double costPerSample = 0.0;

  double cost = 0.0;

  WelfordAggregate() = default;

  WelfordAggregate(SampleCounter ctr, CostInSeconds C, QuantityOfInterest Q,
                   DeltaQuantityOfInterest Y, double cost, double costPerSample) :
      ctr(ctr), C(std::move(C)), Q(std::move(Q)), Y(std::move(Y)),
      cost(cost), costPerSample(costPerSample) {}

  void UpdateTotal();

  void UpdateErrors();

  void UpdateParallel(int commSplit);

  void UpdateCost(double newRoundCost);

  double Value() const { return Q.GetMean(); };

  void UpdateOnComm(const SampleSolution &fSol, const SampleSolution &cSol);

  void UpdateSVarYAndCostPerSample(double sVarY, double expectedCostPerSample);

  friend std::ostream &operator<<(std::ostream &s, const WelfordAggregate &aggregate);
};

#endif //WELFORDAGGREGATE_HPP
