#include "WelfordAggregate.hpp"
#include "Transfers.hpp"
#include "Plotting.hpp"

void Errors::EstimateErrors(const WelfordAggregate &aggregate) {
  disc = EstimateNumeric(aggregate);
  disc2 = disc * disc;
  input = EstimateStochastic(aggregate);
  mse = input + disc2;
  rmse = sqrt(mse);
}

double Errors::EstimateNumeric(const WelfordAggregate &aggregate) {
  if (aggregate.Q.GetMean() == aggregate.Y.GetMean()) return 0.0;
  return std::abs(aggregate.Y.GetMean());
}

double Errors::EstimateStochastic(const WelfordAggregate &aggregate) {
  return aggregate.Q.GetSVar() / double(aggregate.ctr.M);
}

template<typename T>
WelfordDelta<T>::WelfordDelta(const T &newValue) {
  one = newValue;
  ComputeDeltas();
}

template<typename T>
WelfordDelta<T>::WelfordDelta(const T &newValue, const T &oldMean) {
  one = newValue;
  one -= oldMean;
  ComputeDeltas();
}

template<>
WelfordDelta<Vector>::WelfordDelta(const Vector &newValue, const Vector &oldMean) :
    one(newValue), two(0.0, newValue), three(0.0, newValue), four(0.0, newValue) {
  one -= oldMean;
  ComputeDeltas();
}

template<>
WelfordDelta<Vector>::WelfordDelta(const Vector &newValue) :
    one(newValue), two(0.0, newValue), three(0.0, newValue), four(0.0, newValue) {
  ComputeDeltas();
}

template<typename T>
void WelfordDelta<T>::ComputeDeltas() {
  // sOMe FuNny MatHs
  two = one;
  two *= one;
  three = one;
  three *= two;
  four = two;
  four *= two;
}

template<typename T>
WelfordData<T>::WelfordData(
    T mean, T sVar, T skew, T kurt, T powSum2, T powSum3, T powSum4, double weight, double weight2)
    : mean(mean), sVar(sVar), skew(skew), kurt(kurt),
      powSum2(powSum2), powSum3(powSum3), powSum4(powSum4),
      weight(weight), weight2(weight2) {}

template<>
WelfordData<double>::WelfordData(
    double mean, double sVar, double skew, double kurt,
    double powSum2, double powSum3, double powSum4,
    double weight, double weight2)
    : mean(mean), sVar(sVar), kurt(kurt), skew(skew),
      powSum2(powSum2), powSum3(powSum3), powSum4(powSum4),
      weight(weight), weight2(weight2) {}

template<typename T>
WelfordData<T>::WelfordData(const T &newValue)  :
    mean(newValue), sVar(newValue), skew(newValue), kurt(newValue),
    powSum2(newValue), powSum3(newValue), powSum4(newValue), weight(0.0),
    weight2(0.0) {
  mean = 0.0;
  sVar = 0.0;
  skew = 0.0;
  kurt = 0.0;
  powSum2 = 0.0;
  powSum3 = 0.0;
  powSum4 = 0.0;
}

template<>
WelfordData<double>::WelfordData(const WelfordData<double> &welfordData)  :
    mean(welfordData.mean), sVar(welfordData.sVar), skew(welfordData.skew), kurt(welfordData.kurt),
    powSum2(welfordData.powSum2), powSum3(welfordData.powSum3), powSum4(welfordData.powSum4),
    weight(welfordData.weight), weight2(welfordData.weight2) {}

template<>
WelfordData<double>::WelfordData(double value, const WelfordData<double> &welfordData)  :
    mean(value), sVar(value), skew(value), kurt(value),
    powSum2(value), powSum3(value), powSum4(value),
    weight(0.0), weight2(0.0) {}

template<>
WelfordData<Vector>::WelfordData(double value, const WelfordData<Vector> &welfordData)  :
    mean(value, welfordData.mean), sVar(value, welfordData.sVar),
    skew(value, welfordData.skew), kurt(value, welfordData.kurt),
    powSum2(value, welfordData.powSum2), powSum3(value, welfordData.powSum3),
    powSum4(value, welfordData.powSum4),
    weight(0.0), weight2(0.0) {}

template<>
WelfordData<Vector>::WelfordData(double value, const WelfordData<Vector> &welfordData, int commSplit)  :
      mean(Vector(value, welfordData.mean.GetSharedDisc(),
                  welfordData.mean.Level().WithCommSplit(commSplit))),
      sVar(Vector(value, welfordData.sVar.GetSharedDisc(),
                  welfordData.sVar.Level().WithCommSplit(commSplit))),
      skew(Vector(value, welfordData.skew.GetSharedDisc(),
                  welfordData.skew.Level().WithCommSplit(commSplit))),
      kurt(Vector(value, welfordData.kurt.GetSharedDisc(),
                  welfordData.kurt.Level().WithCommSplit(commSplit))),
      powSum2(Vector(value, welfordData.powSum2.GetSharedDisc(),
                     welfordData.powSum2.Level().WithCommSplit(commSplit))),
      powSum3(Vector(value, welfordData.powSum3.GetSharedDisc(),
                     welfordData.powSum3.Level().WithCommSplit(commSplit))),
      powSum4(Vector(value, welfordData.powSum4.GetSharedDisc(),
                     welfordData.powSum4.Level().WithCommSplit(commSplit))),
      weight(welfordData.weight), weight2(welfordData.weight2) {}

template<>
double WelfordData<double>::computeSqrt(const double &value) const {
  return std::sqrt(value);
}

template<>
Vector WelfordData<Vector>::computeSqrt(const Vector &value) const {
  return ComponentSqrt(value);
}


template<typename T>
void Aggregate<T>::UpdateTotal() {
  if (total == nullptr) {
    total = std::make_shared<WelfordData<T>>(0.0, *para);
  }

  total->weight += para->weight;
  total->weight2 += para->weight2;

  WelfordDelta<T> delta(para->mean, total->mean);

  total->computeMean(total->weight - para->weight, para->weight, delta);
  total->computePowSum2(total->weight - para->weight, para->weight, *para, delta);
  total->computePowSum3(total->weight - para->weight, para->weight, *para, delta);
  total->computePowSum4(total->weight - para->weight, para->weight, *para, delta);

  total->computeSVar();
  total->computeSkew();
  total->computeKurt();

  para.reset();
}

template<typename T>
T Aggregate<T>::SumOnCommSplit(const T &value, int commSplit) {
  return PPM->SumOnCommSplit(value, commSplit);
}

template<>
Vector Aggregate<Vector>::SumOnCommSplit(const Vector &value, int commSplit) {
  return SumVectorOnCommSplit(value, commSplit, value.CommSplit());
}

template<typename T>
void Aggregate<T>::UpdateParallel(int maxCommSplit) {
  para = std::make_shared<WelfordData<T>>(0.0, *comm);

  para->weight = PPM->SumAcrossComm(comm->weight, maxCommSplit);
  para->weight2 = PPM->SumAcrossComm(comm->weight2, maxCommSplit);

  if (PPM->Master(maxCommSplit)) {
    para->mean = comm->mean;
    para->powSum2 = comm->powSum2;
    para->powSum3 = comm->powSum3;
    para->powSum4 = comm->powSum4;
  } else {
    para->ResetPowerSums();
  }

  for (int i = maxCommSplit - 1; i >= 0; i--) {
    double W_lr = PPM->SumOnCommSplit(comm->weight, i);
    if (PPM->Master(i)) para->mean *= -1.0;
    WelfordDelta<T> delta(SumOnCommSplit(para->mean, i));
    if (PPM->Master(i)) para->mean *= -1.0;
    double W = PPM->SumOnCommSplit(comm->weight, i + 1);
    para->powSum2 = SumOnCommSplit(para->powSum2, i);
    para->powSum3 = SumOnCommSplit(para->powSum3, i);
    para->powSum4 = SumOnCommSplit(para->powSum4, i);
    if (PPM->Master(i)) {
      para->computeMean((W_lr - W), W, delta);
      para->computePowSum2(W_lr - W, W, delta);
      para->computePowSum3(W_lr - W, W, delta);
      para->computePowSum4(W_lr - W, W, delta);
    } else {
      para->ResetPowerSums();
    }
  }

  if constexpr (std::is_same<T, double>::value) {
    PPM->Broadcast(para->mean);
    PPM->Broadcast(para->powSum2);
    PPM->Broadcast(para->powSum3);
    PPM->Broadcast(para->powSum4);
  } else {
    auto tempPara = std::make_shared<WelfordData<T>>(0.0, *para, 0);
    tempPara->mean = SumVectorAcrossCommSplit(para->mean, maxCommSplit, 0);
    tempPara->powSum2 = SumVectorAcrossCommSplit(para->powSum2, maxCommSplit, 0);
    tempPara->powSum3 = SumVectorAcrossCommSplit(para->powSum3, maxCommSplit, 0);
    tempPara->powSum4 = SumVectorAcrossCommSplit(para->powSum4, maxCommSplit, 0);
    para.reset();
    para = std::make_shared<WelfordData<T>>(*tempPara);
  }

  para->computeSVar();
  para->computeSkew();
  para->computeKurt();

  comm.reset();
}

template<typename T>
void Aggregate<T>::UpdateTemplateOnComm(const T &newValue, double newWeight) {
  bool computeMoments = false;
  if (comm == nullptr) {
    comm = std::make_shared<WelfordData<T>>(newValue);
  } else {
    computeMoments = true;
  }

  comm->weight += newWeight;
  comm->weight2 += newWeight * newWeight;

  WelfordDelta<T> delta(newValue, comm->mean);
  comm->computeMean(comm->weight - newWeight, newWeight, delta);
  comm->computePowSum2(comm->weight - newWeight, newWeight, delta);
  comm->computePowSum3(comm->weight - newWeight, newWeight, delta);
  comm->computePowSum4(comm->weight - newWeight, newWeight, delta);

  if (computeMoments) {
    comm->computeSVar();
    comm->computeSkew();
    comm->computeKurt();
  }
}

void WelfordAggregate::UpdateCost(double newRoundCost) {
  cost += PPM->Max(newRoundCost, 0);
  costPerSample = cost / ctr.M;
}

void WelfordAggregate::UpdateSVarYAndCostPerSample(double sVarY, double expectedCostPerSample) {
  Y = DeltaQuantityOfInterest(WelfordData<double>(0.0, sVarY, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
  costPerSample = expectedCostPerSample;
}

void WelfordAggregate::UpdateOnComm(const SampleSolution &fSol, const SampleSolution &cSol) {
  ctr.Mcomm++;

  C.UpdateOnComm(fSol, cSol);
  Q.UpdateOnComm(fSol, cSol);
  Y.UpdateOnComm(fSol, cSol);
#ifdef AGGREGATE_FOR_SOLUTION
  U.UpdateOnComm(fSol, cSol);
  V.UpdateOnComm(fSol, cSol);
#endif
}

void WelfordAggregate::UpdateParallel(int commSplit) {
  ctr.Mpara = PPM->SumAcrossComm(ctr.Mcomm, commSplit);

  C.UpdateParallel(commSplit);
  Q.UpdateParallel(commSplit);
  Y.UpdateParallel(commSplit);
#ifdef AGGREGATE_FOR_SOLUTION
  U.UpdateParallel(commSplit);
  V.UpdateParallel(commSplit);
#endif

  ctr.Mcomm = 0;
}

void WelfordAggregate::UpdateTotal() {
  ctr.M += ctr.Mpara;

  C.UpdateTotal();
  Q.UpdateTotal();
  Y.UpdateTotal();
#ifdef AGGREGATE_FOR_SOLUTION
  U.UpdateTotal();
  V.UpdateTotal();
#endif

  ctr.Mpara = 0;
}

std::ostream &operator<<(std::ostream &s, const WelfordAggregate &aggregate) {
  s << aggregate.ctr << endl;
  if (aggregate.Q.TotalInitialized())
    s << "QoI  total: " << aggregate.Q.GetTotal();
  if (aggregate.Y.TotalInitialized())
    s << "dQoI total: " << aggregate.Y.GetTotal();
  if (aggregate.C.TotalInitialized())
    s << "cost total: " << aggregate.C.GetTotal();
  if (aggregate.Q.ParaInitialized())
    s << "QoI   para: " << aggregate.Q.GetPara();
  if (aggregate.Y.ParaInitialized())
    s << "dQoI  para: " << aggregate.Y.GetPara();
  if (aggregate.C.ParaInitialized())
    s << "cost  para: " << aggregate.C.GetPara();
  if (aggregate.Q.CommInitialized())
    s << "QoI   comm: " << aggregate.Q.GetComm();
  if (aggregate.Y.CommInitialized())
    s << "dQoI  comm: " << aggregate.Y.GetComm();
  if (aggregate.C.CommInitialized())
    s << "cost  comm: " << aggregate.C.GetComm();
  return s;
}