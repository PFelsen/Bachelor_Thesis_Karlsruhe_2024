#include "Parallel.hpp"
#include "ExchangeBuffer.hpp"
#include "MemoryLogger.hpp"

#include <cmath>
#include <iostream>


#ifdef BUILD_IA

#include "IACInterval.hpp"


#endif

ParallelProgrammingModel *ParallelProgrammingModel::ppm = nullptr;

ParallelProgrammingModel::ParallelProgrammingModel(int *argc, char ***argv) {
  MPI_Init(argc, argv);

#ifdef BUILD_IA
  MPI_Op_create((MPI_User_function *)mpi_sum_iainterval, 1, &MPI_SUM_IAINTERVAL);
  MPI_Op_create((MPI_User_function *)mpi_sum_iacinterval, 1, &MPI_SUM_IACINTERVAL);
  MPI_Op_create((MPI_User_function *)mpi_hull_iainterval, 1, &MPI_HULL_IAINTERVAL);
  MPI_Op_create((MPI_User_function *)mpi_intersect_iainterval, 1, &MPI_INTERSECT_IAINTERVAL);
#endif

  comms.push_back(new WorldCommunicator());
}

ParallelProgrammingModel::~ParallelProgrammingModel() {
  if (MasterLogging::isInitialized()) mout.CommunicateWarnings();

#ifdef BUILD_IA
  MPI_Op_free(&MPI_SUM_IAINTERVAL);
  MPI_Op_free(&MPI_SUM_IACINTERVAL);
  MPI_Op_free(&MPI_HULL_IAINTERVAL);
  MPI_Op_free(&MPI_INTERSECT_IAINTERVAL);
#endif

  ClearCommunicators();
  delete comms[0];
  comms.clear();
  MPI_Finalize();
}

void ParallelProgrammingModel::Initialize(int argc, char **argv) {
  if (IsInitialized()) THROW("PPM already initialized!")
  ppm = new ParallelProgrammingModel(&argc, &argv);
  PPM->FullSplit();
}

bool ParallelProgrammingModel::IsInitialized() { return ppm != nullptr; }

ParallelProgrammingModel *ParallelProgrammingModel::Instance() {
  if (!ppm) THROW("PPM not initialized!")
  return ppm;
}

void ParallelProgrammingModel::Close() {
  if (ppm) delete ppm;
  ppm = nullptr;
}

void ParallelProgrammingModel::Abort() { MPI_Abort(MPI_COMM_NULL, 1); }

Communicator *ParallelProgrammingModel::GetComm(int commSplit) const {
  if (commSplit == -1) return comms.back();
  if (commSplit > MaxCommSplit()) {
    // TODO: Specify behaviour
    return comms.back();
    THROW("Communicator split " + std::to_string(commSplit) + " does not exist")
  }
  return comms.at(commSplit);
}

Communicator *ParallelProgrammingModel::Copy(int commSplit) const {
  if (commSplit == -1) return comms.back()->Copy();
  if (commSplit > MaxCommSplit()) {
    THROW("Communicator split" + std::to_string(commSplit) + " does not exist")
  }
  return comms.at(commSplit)->Copy();
}

int ParallelProgrammingModel::Proc(int commSplit) const { return GetComm(commSplit)->Proc(); }

int ParallelProgrammingModel::Size(int commSplit) const { return GetComm(commSplit)->Size(); }

bool ParallelProgrammingModel::Master(int commSplit) const { return GetComm(commSplit)->Master(); }

int ParallelProgrammingModel::Color(int commSplit) const { return GetComm(commSplit)->Color(); }

bool ParallelProgrammingModel::IsParallel(int commSplit) const { return Size(commSplit) > 1; }

void ParallelProgrammingModel::PrintInfo(int commSplit) const {
  return GetComm(commSplit)->PrintInfo();
}

void ParallelProgrammingModel::Communicate(ExchangeBuffer &exBuffer, int commSplit) {
  GetComm(commSplit)->Communicate(exBuffer);
}

bool ParallelProgrammingModel::And(bool b, int commSplit) const {
  return GetComm(commSplit)->And(b);
}

bool ParallelProgrammingModel::Or(bool b, int commSplit) const { return GetComm(commSplit)->Or(b); }

bool ParallelProgrammingModel::SplitCommunicators() {
  PPM->Barrier(0);
  bool fullSplitReached = true;
  auto comm = comms.back();
  for (int p = 0; p < PPM->Size(0); p++) {
    if (p == PPM->Proc(0) && comm->Size() != 1) {
      fullSplitReached = false;
      break;
    }
  }
  if (!fullSplitReached) comms.push_back(comm->Split());
  return fullSplitReached;
}

bool ParallelProgrammingModel::SplitCommunicators(int commSplit) {
  if (commSplit == -1) {
    FullSplit();
    return true;
  }
  bool fullSplitReached = false;
  for (int i = 0; i < commSplit; i++)
    fullSplitReached = SplitCommunicators();
  return fullSplitReached;
}

void ParallelProgrammingModel::FullSplit() {
  bool fullSplitReached = false;
  while (!fullSplitReached)
    fullSplitReached = SplitCommunicators();
}

int ParallelProgrammingModel::CommsSize() const { return comms.size(); }

int ParallelProgrammingModel::MaxCommSplit() const {
  if (CommsSize() != 0) return CommsSize() - 1;
  else return 0;
}

void ParallelProgrammingModel::ClearCommunicators() {
  while (comms.size() > 1) {
    delete comms[comms.size() - 1];
    comms.pop_back();
  }
}

void ParallelProgrammingModel::Barrier(int commSplit) { GetComm(commSplit)->Barrier(); }

#ifdef BUILD_IA
void ParallelProgrammingModel::Hull(IAInterval *a, size_t n, int commSplit) {

  GetComm(commSplit)->Hull(a, n);
}

void ParallelProgrammingModel::Intersect(IAInterval *a, size_t n, int commSplit) {
  GetComm(commSplit)->Intersect(a, n);
}

IAInterval ParallelProgrammingModel::Hull(IAInterval a) {
  Hull(&a, 1);
  return a;
}

IAInterval ParallelProgrammingModel::Intersect(IAInterval a) {
  Intersect(&a, 1);
  return a;
}

#endif


// Credits to: https://gist.github.com/thirdwing/da4621eb163a886a03c5
#include <fstream>
#include <unistd.h>

void ParallelProgrammingModel::Synchronize(const char *c) {
  std::cout.flush();
  PPM->Barrier(0);
  if (c == 0) return;
  std::cout << GetComm(-1)->Proc() << "|";
  PPM->Barrier(0);
  std::cout.flush();
  if (PPM->Master(0)) std::cout << c << endl;
  std::cout.flush();
  PPM->Barrier(0);
  std::cout.flush();
}

/*
 * =========================
 * Deprecated
 * =========================
 */

bool ParallelProgrammingModel::Boolean(bool b) { return GetComm(0)->Or(b); }

ParallelProgrammingModel *ParallelProgrammingModel::instance() {
  if (!ppm) THROW("PPM not initialized!")
  return ppm;
}

bool ParallelProgrammingModel::isInitialized() { return ppm != nullptr; }

void ParallelProgrammingModel::BroadcastDouble(double d) {
  Assert(Master(0));
  BCastOnCommSplit(d, 0);
}

double ParallelProgrammingModel::TotalMemoryUsage(bool perPro) {
  return MemoryLogger::TotalMemoryUsage(perPro);
}

double ParallelProgrammingModel::BroadcastDouble() {
  Assert(!Master(0));
  double a;
  BCastOnCommSplit(a, 0);
  return a;
}

void ParallelProgrammingModel::BroadcastInt(int i) {
  Assert(Master(0));
  BCastOnCommSplit(i, 0);
}

int ParallelProgrammingModel::BroadcastInt() {
  Assert(!Master(0));
  int i;
  BCastOnCommSplit(i, 0);
  return i;
}

bool ParallelProgrammingModel::master() const { return GetComm(0)->Master(); }

int ParallelProgrammingModel::proc() const { return GetComm(0)->Proc(); }

int ParallelProgrammingModel::size() const { return GetComm(0)->Size(); }

void ParallelProgrammingModel::close() { Close(); }

void ParallelProgrammingModel::abort() { MPI_Abort(MPI_COMM_NULL, 1); }

void ParallelProgrammingModel::initialize(int *argc, char **argv) { Initialize(*argc, argv); }

int ParallelProgrammingModel::NumberOfCommunicators(int commSplit) {
  return int(pow(2, commSplit));
}
