#include "Assertion.hpp"

#include <algorithm>
#include <iostream>

#include "Config.hpp"
#include "Parallel.hpp"

std::vector<std::string> Assertion::previousWarnings = {};

void Assertion::assertion(bool assertion, const std::string &message, const std::string &file,
                          size_t line) {
  if (assertion) return;
  std::string msg = "\n\033[1;35mAssert: " + message + "\033[0m\n";
  if (ParallelProgrammingModel::IsInitialized())
    msg += "        on proc " + std::to_string(PPM->Proc(0)) + "\n";
  const auto lineString = std::to_string(line);
  msg += "        in " + file + ":" + lineString + "\n";
  msg += "        on line " + lineString + "\n\n";
  std::cerr << msg;
  std::cerr << std::flush;
  assert(0);
}

void Assertion::assertion(bool b, const char *s, int n, const char *f, int L) {
  if (b) return;
  std::string msg = "\n\033[1;35mAssert: " + std::string(s) + std::to_string(n) + "\033[0m\n";
  if (ParallelProgrammingModel::IsInitialized())
    msg += "        on proc " + std::to_string(PPM->Proc(0)) + "\n";
  msg += "        in " + std::string(f) + ":" + std::to_string(L) + "\n";
  msg += "        on line " + std::to_string(L) + "\n\n";
  std::cerr << msg;
  std::cerr << std::flush;
  assert(0);
}

void Assertion::exit(const char *s, const char *f, int L) { exit(std::string(s), f, L); }

void Assertion::exit(std::string s, const char *f, int L) {
  std::string msg = "\n\033[1;31mError: " + s + "\033[0m\n";
  if (ParallelProgrammingModel::IsInitialized())
    msg += "       on proc " + std::to_string(PPM->Proc(0)) + "\n";
  msg += "       in " + std::string(f) + ":" + std::to_string(L) + "\n";
  msg += "       on line " + std::to_string(L) + "\n\n";
  std::cerr << msg;
  std::cerr << std::flush;
}

void Assertion::warning(const char *s, const char *f, int L) { warning(std::string(s), f, L); }

void Assertion::warning(std::string s, const char *f, int L) {
  if (!PPM->Master(0)) return;

  if (std::find(previousWarnings.begin(), previousWarnings.end(), s) != previousWarnings.end())
    return;
  previousWarnings.push_back(s);

  std::string msg = "\n\033[1;33mWarning: " + s + "\033[0m\n" + "         in " + std::string(f)
                    + ":" + std::to_string(L) + "\n" + "         on line " + std::to_string(L)
                    + "\n\n";
  std::cout << msg;
  std::cout << std::flush;
  mout.AddWarningMsg("Warning: " + s + "\n" + "         in " + std::string(f) + ":"
                     + std::to_string(L) + "\n" + "         on line " + std::to_string(L) + "\n\n");
}

void Assertion::warningOnProc(int proc, const char *s, const char *f, int L) {
  warningOnProc(proc, std::string(s), f, L);
}

void Assertion::warningOnProc(int proc, std::string s, const char *f, int L) {
  if (PPM->Proc(0) != proc) return;

  if (std::find(previousWarnings.begin(), previousWarnings.end(), s) != previousWarnings.end())
    return;
  previousWarnings.push_back(s);

  std::string msg = "\n\033[1;33mWarning: " + s + "\033[0m\n";
  if (ParallelProgrammingModel::IsInitialized())
    msg += "         on proc " + std::to_string(PPM->Proc(0)) + "\n";
  msg += "         in " + std::string(f) + ":" + std::to_string(L) + "\n";
  msg += "         on line " + std::to_string(L) + "\n\n";
  std::cout << msg;
  std::cout << std::flush;

  msg = "Warning: " + s + "\n";
  if (ParallelProgrammingModel::IsInitialized())
    msg += "         on proc " + std::to_string(PPM->Proc(0)) + "\n";
  msg += "         in " + std::string(f) + ":" + std::to_string(L) + "\n";
  msg += "         on line " + std::to_string(L) + "\n\n";
  mout.AddWarningMsg(msg);
}

void Assertion::warningOnProc(const char *s, const char *f, int L) {
  warningOnProc(std::string(s), f, L);
}

void Assertion::warningOnProc(std::string s, const char *f, int L) {
  if (std::find(previousWarnings.begin(), previousWarnings.end(), s) != previousWarnings.end())
    return;
  previousWarnings.push_back(s);

  std::string msg = "\n\033[1;33mWarning: " + s + "\033[0m\n";
  if (ParallelProgrammingModel::IsInitialized())
    msg += "         on proc " + std::to_string(PPM->Proc(0)) + "\n";
  msg += "         in " + std::string(f) + ":" + std::to_string(L) + "\n";
  msg += "         on line " + std::to_string(L) + "\n\n";
  std::cout << msg;
  std::cout << std::flush;

  msg = "Warning: " + s + "\n";
  if (ParallelProgrammingModel::IsInitialized())
    msg += "         on proc " + std::to_string(PPM->Proc(0)) + "\n";
  msg += "         in " + std::string(f) + ":" + std::to_string(L) + "\n";
  msg += "         on line " + std::to_string(L) + "\n\n";
  mout.AddWarningMsg(msg);
}

void Assertion::synchronizeErrors(bool failed, std::string s, const char *f, int L) {
  std::string msg = "";
  if (failed) {
    exit(s.c_str(), f, L);
    msg += "Error: " + s + "\n";
    if (ParallelProgrammingModel::IsInitialized())
      msg += "       on proc " + std::to_string(PPM->Proc(0)) + "\n";
    msg += "       in " + std::string(f) + ":" + std::to_string(L) + "\n";
    msg += "       on line " + std::to_string(L) + "\n\n";
  }
  failed = PPM->Or(failed, 0);
  if (failed) {
    const char *c_msg = msg.c_str();
    ExchangeBuffer exBuffer(0);
    for (int i = 0; i < msg.length(); ++i)
      exBuffer.Send(0) << c_msg[i];
    exBuffer.Communicate();
    if (PPM->Proc(0) == 0) {
      for (short p = 0; p < PPM->Size(0); ++p) {
        std::string msg_p;
        while (exBuffer.Receive(p).size() < exBuffer.ReceiveSize(p)) {
          char c;
          exBuffer.Receive(p) >> c;
          msg_p += c;
        }
        mout.AddErrorMsg(msg_p);
      }
    }
    std::exit(1);
  }
}

void Assertion::Error(std::string s, const char *f, int L) {
  std::string msg = "";
  exit(s.c_str(), f, L);
  msg += "Error: " + s + "\n";
  if (ParallelProgrammingModel::IsInitialized()) msg += "       on all procs\n";
  msg += "       in " + std::string(f) + ":" + std::to_string(L) + "\n";
  msg += "       on line " + std::to_string(L) + "\n\n";

  if (PPM->Proc(0) == 0) { mout.AddErrorMsg(msg); }
}
