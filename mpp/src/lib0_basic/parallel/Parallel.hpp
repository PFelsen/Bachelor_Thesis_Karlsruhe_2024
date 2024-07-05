#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "Assertion.hpp"
#include "Communicator.hpp"

class ParallelProgrammingModel {
private:
  static ParallelProgrammingModel *ppm;
public:
  static void Initialize(int argc, char **argv);

  static bool IsInitialized();

  static ParallelProgrammingModel *Instance();

  static void Close();

  static void Abort();

  static int NumberOfCommunicators(int commSplit = 0);
private:
  Communicators comms;
public:
  ParallelProgrammingModel(int *argc, char ***argv);

  ~ParallelProgrammingModel();

  int Proc(int commSplit = 0) const;

  int Size(int commSplit = 0) const;

  bool Master(int commSplit = 0) const;

  int Color(int commSplit = 0) const;

  bool IsParallel(int commSplit) const;

  void PrintInfo(int commSplit = 0) const;

  Communicator *GetComm(int commSplit = 0) const;

  Communicator *Copy(int commSplit = 0) const;

  // TODO: Remove size issue such that the following holds (NOT true yet!!!):
  /// Note that size is only the size of the array (NOT multiplied with sizeof(T))!
  template<typename T>
  void Broadcast(T *data, size_t size, int commSplit = 0, int q = 0) {
    size /= sizeof(T);
    GetComm(commSplit)->Broadcast(data, size, q);
  }

  template<typename T>
  void BCastOnCommSplit(T &t, int commSplit) {
    GetComm(commSplit)->Broadcast(t);
  }

  template<typename T>
  void Broadcast(const T &t) {
    GetComm(0)->Broadcast(t);
  }

  template<typename T>
  void BCastOnCommSplit(const T &t, int commSplit) {
    GetComm(commSplit)->Broadcast(t);
  }

  template<typename T>
  void Send(T *data, size_t size, int destination, int commSplit = 0) const {
    GetComm(commSplit)->Send(data, size, destination);
  };

  template<typename T>
  void Receive(T *data, size_t size, int source, int commSplit = 0) const {
    GetComm(commSplit)->Receive(data, size, source);
  };

  // Tested in Exchange Buffer
  void Communicate(ExchangeBuffer &exBuffer, int commSplit = 0);

  template<typename T>
  void Allgather(T *sendData, size_t sizeSend, T *receiveData, int commSplit = 0) const {
    GetComm(commSplit)->Allgather(sendData, sizeSend, receiveData);
  };

  template<typename T>
  void Allgatherv(T *sendData, size_t sizeSend, T *receiveData, int *sizeReceive, int *displs,
                  int commSplit) const {
    GetComm(commSplit)->Allgatherv(sendData, sizeSend, receiveData, sizeReceive, displs);
  };

  template<typename T>
  void Sum(T *a, size_t n, int commSplit = 0) const {
    GetComm(commSplit)->Sum(a, n);
  };

  template<typename T>
  T SumOnCommSplit(T a, int commSplit) const {
    Sum(&a, 1, commSplit);
    return a;
  }

  template<typename T>
  void Sum(std::vector<T> &v) const {
    Sum(v.data(), std::size(v));
  }

  template<typename T>
  void SumOnCommSplit(std::vector<T> &v, int commSplit) const {
    size_t n = std::size(v);
    T *a = new T[n];
    for (int i = 0; i < n; ++i)
      a[i] = v[i];
    Sum(a, n, commSplit);
    for (int i = 0; i < n; ++i)
      v[i] = a[i];
    delete[] a;
  }

  template<typename T>
  T SumAcrossComm(T a, int commSplit) const {
    return (SumOnCommSplit(a, 0) / Size(commSplit));
  }

  template<typename T>
  T Min(T a, int commSplit = 0) const {
    return GetComm(commSplit)->Min(a);
  }

  template<typename T>
  T Max(T a, int commSplit = 0) const {
    return GetComm(commSplit)->Max(a);
  }

  template<typename T>
  T ExclusiveScanSum(T a, int commSplit = 0) const {
    return GetComm(commSplit)->ExclusiveScanSum(a);
  }

  bool And(bool b, int commSplit = 0) const;

  bool Or(bool b, int commSplit = 0) const;

#ifdef BUILD_IA

  void Hull(IAInterval *a, size_t n, int commSplit = 0);

  void Intersect(IAInterval *a, size_t n, int commSplit = 0);

  IAInterval Hull(IAInterval a);

  IAInterval Intersect(IAInterval a);

#endif

  bool SplitCommunicators();

  bool SplitCommunicators(int commSplit);

  void FullSplit();

  int CommsSize() const;

  int MaxCommSplit() const;

  /// Deletes a communicators except world
  void ClearCommunicators();

  void Barrier(int commSplit = 0);

  void Synchronize(const char *c = 0);

  double TotalMemoryUsage(bool perPro = false);

  void BroadcastDouble(double d);

  double BroadcastDouble();

  void BroadcastInt(int i);

  int BroadcastInt();

  /*
   * Functions below should not be used. Instead the capitalized version of the function with
   * CommSplit should be used
   */

  template<typename T>
  void Broadcast(T &t) {
    GetComm(0)->Broadcast(t);
  }

  template<typename T>
  T Sum(T a) const {
    Sum(&a, 1);
    return a;
  }

  bool Boolean(bool b);

  static bool isInitialized();

  static void initialize(int *argc, char **argv);

  static ParallelProgrammingModel *instance();

  int proc() const;

  int size() const;

  bool master() const;

  static void abort();

  static void close();
};

#define PPM (ParallelProgrammingModel::Instance())

#endif // of #ifndef _PARALLEL_H_
