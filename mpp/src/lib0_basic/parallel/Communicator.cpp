#include "Communicator.hpp"

#include <iostream>

#ifdef BUILD_IA

#include "IACInterval.hpp"

#endif

Communicator::Communicator(MPI_Comm comm, int commSplit, int color) :
    P(0), p(0), mpiComm(comm), commSplit(commSplit), color(color) {
  MPI_Comm_size(comm, &P);
  MPI_Comm_rank(comm, &p);
}

Communicator::Communicator(const Communicator &comm, int color) :
    P(0), p(0), commSplit(-1), color(color) {
  MPI_Comm_split(comm.mpiComm, color, comm.p, &mpiComm);
  MPI_Comm_size(mpiComm, &P);
  MPI_Comm_rank(mpiComm, &p);
}

Communicator::~Communicator() {
  if (mpiComm != MPI_COMM_WORLD) MPI_Comm_free(&mpiComm);
}

bool Communicator::Master() const { return (p == 0); }

int Communicator::Proc() const { return p; }

int Communicator::Size() const { return P; }

int Communicator::Color() const { return color; }

int Communicator::CommSplit() const { return commSplit; }

Communicator *Communicator::Split() const {
  int worldSize;
  int worldRank;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  int newColor = worldRank * pow(2, commSplit + 1) / worldSize;
  MPI_Comm newComm;
  MPI_Comm_split(mpiComm, newColor, p, &newComm);
  return new Communicator(newComm, commSplit + 1, newColor);
}

Communicator *Communicator::Copy() const {
  MPI_Comm newComm;
  MPI_Comm_dup(mpiComm, &newComm);
  return new Communicator(newComm, commSplit, color);
}

#ifdef BUILD_IA

void Communicator::Sum(IAInterval *a, size_t n) const {
  double in[2 * n];
  double inout[2 * n];
  for (int i = 0; i < n; ++i) {
    in[2 * i] = a[i].inf();
    in[2 * i + 1] = a[i].sup();
  }
  MPI_Allreduce(in, inout, 2 * n, MPI_DOUBLE, MPI_SUM_IAINTERVAL, mpiComm);
  for (int i = 0; i < n; ++i)
    a[i] = IAInterval(inout[2 * i], inout[2 * i + 1]);
}

void Communicator::Sum(IACInterval *a, size_t n) const {
  double in[4 * n];
  double inout[4 * n];
  for (int i = 0; i < n; ++i) {
    in[4 * i] = a[i].real().inf();
    in[4 * i + 1] = a[i].real().sup();
    in[4 * i + 2] = a[i].imag().inf();
    in[4 * i + 3] = a[i].imag().sup();
  }
  MPI_Allreduce(in, inout, 4 * n, MPI_DOUBLE, MPI_SUM_IACINTERVAL, mpiComm);
  for (int i = 0; i < n; ++i)
    a[i] = IACInterval(IAInterval(inout[4 * i], inout[4 * i + 1]),
                       IAInterval(inout[4 * i + 2], inout[4 * i + 3]));
}

void Communicator::Hull(IAInterval *a, size_t n) const {
  double in[2 * n];
  double inout[2 * n];
  for (int i = 0; i < n; ++i) {
    in[2 * i] = a[i].inf();
    in[2 * i + 1] = a[i].sup();
  }
  MPI_Allreduce(in, inout, 2 * n, MPI_DOUBLE, MPI_HULL_IAINTERVAL, mpiComm);
  for (int i = 0; i < n; ++i)
    a[i] = IAInterval(inout[2 * i], inout[2 * i + 1]);
}

void Communicator::Intersect(IAInterval *a, size_t n) const {
  double in[2 * n];
  double inout[2 * n];
  for (int i = 0; i < n; ++i) {
    in[2 * i] = a[i].inf();
    in[2 * i + 1] = a[i].sup();
  }
  MPI_Allreduce(in, inout, 2 * n, MPI_DOUBLE, MPI_INTERSECT_IAINTERVAL, mpiComm);
  for (int i = 0; i < n; ++i)
    a[i] = IAInterval(inout[2 * i], inout[2 * i + 1]);
}

#endif

void Communicator::PrintInfo() const {
  std::cout << "\tProc rank: " << p << "\tComm size: " << P << "\tComm color: " << color
            << "\tComm split: " << commSplit << std::endl;
}

void Communicator::Communicate(ExchangeBuffer &exBuffer) {
  if (exBuffer.SendMessages() == 0 && exBuffer.RecvMessages() == 0) { return; }
  const int tag = 27;
  MPI_Request sendRequest[exBuffer.SendMessages()];
  MPI_Request recvRequest[exBuffer.RecvMessages()];
  MPI_Request *r = recvRequest;
  for (int q = 0; q < Size(); ++q)
    for (int k = exBuffer.Messages(q); k < exBuffer.Messages(q + 1); ++k)
      if (exBuffer.MessageDest(k) == Proc())
        MPI_Irecv(exBuffer.Receive(q)(), exBuffer.MessageSize(k), MPI_BYTE, q, tag, mpiComm, r++);
  r = sendRequest;
  for (int k = exBuffer.Messages(Proc()); k < exBuffer.Messages(Proc() + 1); ++k) {
    MPI_Isend(exBuffer.Send(exBuffer.MessageDest(k))(), exBuffer.MessageSize(k), MPI_BYTE,
              exBuffer.MessageDest(k), tag, mpiComm, r++);
  }
  MPI_Status st;
  r = sendRequest;
  for (int k = 0; k < exBuffer.SendMessages(); ++k)
    MPI_Wait(r++, &st);
  r = recvRequest;
  for (int k = 0; k < exBuffer.RecvMessages(); ++k)
    MPI_Wait(r++, &st);
}

bool Communicator::And(bool b) {
  int a = 0;
  if (b) a = 1;
  a = Min(a);
  if (a == 0) return false;
  return true;
}

bool Communicator::Or(bool b) {
  int a = 0;
  if (b) a = 1;
  a = Max(a);
  if (a == 1) return true;
  return false;
}

void Communicator::Barrier() {
  // Blocks the caller until all group members have called it.
  // The call returns at any process only after all group
  // members have entered the call.
  MPI_Barrier(mpiComm);
}
