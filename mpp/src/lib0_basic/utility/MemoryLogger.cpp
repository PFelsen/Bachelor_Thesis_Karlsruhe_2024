#include "MemoryLogger.hpp"

int MemoryLogger::verbose = 1;

double MemoryLogger::initTime = 0.0;

double MemoryLogger::offSetMemory = 0.0;

std::vector<MemoryStamp> MemoryLogger::memoryStamps = {};

std::vector<double> MemoryLogger::GetTimeStampVector() {
  std::vector<double> dates;
  for (const auto &stamp : memoryStamps)
    dates.push_back(stamp.timeSinceTareInSeconds);
  return dates;
}

std::vector<double> MemoryLogger::GetMemoryStampVector() {
  std::vector<double> memories;
  for (const auto &stamp : memoryStamps)
    memories.push_back(stamp.memorySinceTareInMByte);
  return memories;
}

std::vector<std::string> MemoryLogger::GetDescriptionVector() {
  std::vector<std::string> memories;
  for (const auto &stamp : memoryStamps)
    memories.push_back(stamp.description);
  return memories;
}

void MemoryLogger::LogMemory(const std::string &description) {
  double memUsageSinceTare = MemoryUsage(false) - offSetMemory;
  vout(2) << description << " mem=" << memUsageSinceTare << " MB";
  vout(3) << " (since tare)";
  vout(2) << endl;
  memoryStamps.emplace_back(description, memUsageSinceTare, MPI_Wtime() - initTime);
}

void MemoryLogger::Tare() {
  Config::Get("MemoryLoggerVerbose", verbose);
  offSetMemory = MemoryUsage();
  initTime = MPI_Wtime();
  memoryStamps = {};
}

void MemoryLogger::PrintInfo() {
  mout.PrintInfo("Memory", verbose, PrintInfoEntry("Time-stamp", GetTimeStampVector(), 1),
                 PrintInfoEntry("Memory", GetMemoryStampVector(), 1),
                 PrintInfoEntry("Description", GetDescriptionVector(), 1));
}

void processMemoryUsage(double &vm_usage, double &resident_set) {
  vm_usage = 0.0;
  resident_set = 0.0;

  // the two fields we want
  unsigned long vsize;
  long rss;
  {
    std::string ignore;
    std::ifstream ifs("/proc/self/stat", std::ios_base::in);
    ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> ignore >> ignore >> vsize >> rss;
    ifs.close();
  }

  // in case x86-64 is configured to use 2MB pages
  long page_size_b = sysconf(_SC_PAGE_SIZE);
  vm_usage = (double)vsize / (1024.0 * 1024.0);
  resident_set = (double)(rss * page_size_b) / (1024.0 * 1024.0);
}

double MemoryLogger::MemoryUsage(bool perProcess) {
  // returns the physical memory im MB for all processes
  // (excluding shared libraries and swap)
  double vm_usage, resident_set;
  PPM->Synchronize();
  processMemoryUsage(vm_usage, resident_set);

  if (perProcess) {
    ExchangeBuffer exBuffer(0);

    if (!PPM->Master(0)) {
      exBuffer.Send(0) << PPM->Proc(0);
      exBuffer.Send(0) << resident_set;
      exBuffer.Communicate();
    } else {
      exBuffer.Communicate();

      mout << "    proc 0 mem used: " << resident_set << " MB" << endl;
      for (short q = 0; q < PPM->Size(0); ++q) {
        while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
          int p = -1;
          double p_resident_set = -1;

          exBuffer.Receive(q) >> p;
          exBuffer.Receive(q) >> p_resident_set;

          mout << "    proc " << p << " mem used: " << p_resident_set << " MB" << endl;
        }
      }
    }
  }

  resident_set = PPM->SumOnCommSplit(resident_set, 0);
  return resident_set;
}

double MemoryLogger::TotalMemoryUsage(bool perProcess) {
  double vm_usage, resident_set;
  PPM->Synchronize();
  processMemoryUsage(vm_usage, resident_set);

  if (perProcess)
    pout << "Memory on proc " << PPM->Proc(0) << " is: " << resident_set << " MB" << endl;

  resident_set = PPM->SumOnCommSplit(resident_set, 0);
  return resident_set;
}