#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include <type_traits>
#include "ExchangeBuffer.hpp"

#include "DataTypes.hpp"

#ifdef BUILD_IA

#include "IAMpiOperators.hpp"

#endif

#include <complex>
#include <vector>


constexpr size_t MaxBroadcastSize = 2147483647;

// 2^31 - 1

class Communicator {
protected:
  // communicator size
  int P;

  // process rank
  int p;

  // communicator color
  int color;

  // communicators are divided into 2^commSplit partitions
  int commSplit;

  MPI_Comm mpiComm;
public:
  Communicator(MPI_Comm comm, int commSplit, int color);

  Communicator(const Communicator &comm, int color);

  ~Communicator();

  bool Master() const;

  int Proc() const;

  int Size() const;

  int Color() const;

  int CommSplit() const;

  virtual void PrintInfo() const;

  Communicator *Split() const;

  Communicator *Copy() const;

  /// Broadcast for known mpi data types
  /// Note that size is only the size of the array (NOT multiplied with sizeof(T))!
  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, void>::type Broadcast(T *data, size_t size,
                                                                            int q = 0) {
    while (size > MaxBroadcastSize) {
      MPI_Bcast(data, MaxBroadcastSize, mpi_datatype<T>::value(), q, mpiComm);
      MPI_Barrier(mpiComm);
      data += MaxBroadcastSize;
      size -= MaxBroadcastSize;
    }

    if (size > 0) { MPI_Bcast(data, size, mpi_datatype<T>::value(), q, mpiComm); }
    MPI_Barrier(mpiComm);
  }

  /// Broadcast for unknown mpi data types via char array
  /// Note that size is only the size of the array (NOT multiplied with sizeof(T))!
  template<typename T>
  typename std::enable_if<!has_mpi_datatype<T>::value, void>::type Broadcast(T *data, size_t size,
                                                                             int q = 0) {
    char *Data = (char *)data;
    size *= sizeof(T);

    while (size > MaxBroadcastSize) {
      MPI_Bcast(Data, MaxBroadcastSize, MPI_BYTE, q, mpiComm);
      MPI_Barrier(mpiComm);
      Data += MaxBroadcastSize;
      size -= MaxBroadcastSize;
    }

    if (size > 0) { MPI_Bcast(Data, size, MPI_BYTE, q, mpiComm); }
    MPI_Barrier(mpiComm);
  }

  template<typename T>
  void Broadcast(T &t) {
    Broadcast(&t, 1);
  }

  template<typename T>
  void Broadcast(const T &t) {
    Assert(Master());
    T data = t;
    Broadcast(&data, 1);
  }

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, void>::type Send(T *data, size_t size,
                                                                       int destination) const {
    const int tag = 27; // TODO: Why 27?
    MPI_Request request;
    MPI_Isend(data, size, mpi_datatype<T>::value(), destination, tag, mpiComm, &request);
    MPI_Status status;
    MPI_Wait(&request, &status);
  }

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, void>::type Receive(T *data, size_t size,
                                                                          int source) const {
    const int tag = 27; // TODO: Why 27?
    MPI_Request request;
    MPI_Irecv(data, size, mpi_datatype<T>::value(), source, tag, mpiComm, &request);
    MPI_Status status;
    MPI_Wait(&request, &status);
  }

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, void>::type
  Allgather(T *sendData, size_t sizeSend, T *receiveData) const {
    MPI_Allgather(sendData, sizeSend, mpi_datatype<T>::value(), receiveData, sizeSend,
                  mpi_datatype<T>::value(), mpiComm);
  }

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, void>::type
  Allgatherv(T *sendData, size_t sizeSend, T *receiveData, int *sizeReceive, int *displs) const {
    MPI_Allgatherv(sendData, sizeSend, mpi_datatype<T>::value(), receiveData, sizeReceive, displs,
                   mpi_datatype<T>::value(), mpiComm);
  }

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, void>::type Sum(T *a, size_t n) const {
    if (n > 0) {
      T b[n];
      MPI_Allreduce(a, b, n, mpi_datatype<T>::value(), MPI_SUM, mpiComm);
      memcpy(a, b, sizeof(T) * n);
    }
  }

  template<template<typename, int, int> class P, typename T, int sDim, int tDim>
  typename std::enable_if<has_mpi_datatype<T>::value, void>::type Sum(P<T, sDim, tDim> *p,
                                                                      size_t n) const {
    T a[n * (sDim + tDim)];

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < sDim + tDim; ++j)
        a[n * i + j] = p[i][j];
    Sum(a, n * (sDim + tDim));
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < sDim + tDim; ++j)
        p[i][j] = a[n * i + j];
  }

#ifdef BUILD_IA
  void Sum(IAInterval *a, size_t n) const;

  void Sum(IACInterval *a, size_t n) const;

  void Hull(IAInterval *a, size_t n) const;

  void Intersect(IAInterval *a, size_t n) const;
#endif

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, T>::type Min(T a) {
    T b;
    MPI_Allreduce(&a, &b, 1, mpi_datatype<T>::value(), MPI_MIN, mpiComm);
    return b;
  };

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, T>::type Max(T a) {
    T b;
    MPI_Allreduce(&a, &b, 1, mpi_datatype<T>::value(), MPI_MAX, mpiComm);
    return b;
  };

  void Communicate(ExchangeBuffer &exBuffer);

  bool And(bool b);

  bool Or(bool b);

  void Barrier();

  template<typename T>
  typename std::enable_if<has_mpi_datatype<T>::value, T>::type ExclusiveScanSum(T a) {
    T b;
    MPI_Exscan(&a, &b, 1, mpi_datatype<T>::value(), MPI_SUM, mpiComm);
    return b;
  }
};

struct WorldCommunicator : public Communicator {
  WorldCommunicator() : Communicator(MPI_COMM_WORLD, 0, 0){};
};

typedef std::vector<Communicator *> Communicators;

#endif // COMMUNICATOR_HPP
