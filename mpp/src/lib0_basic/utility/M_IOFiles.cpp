#include <cstring>
#include <iostream>

#include "M_IOFiles.hpp"
#include "Parallel.hpp"
#include "ctools.hpp"


using namespace std;

const string MError("M_ERROR cannot open file ");

char namebuffer[64];

// char NameBuffer[64];

M_ifstream::M_ifstream(const char *name, bool test) : std::ifstream(name) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  if (test)
    if (!*this) {
      std::cerr << "M_ERROR cannot open file " << name << endl;
      exit(1);
    }
}

M_ofstream::M_ofstream(const char *name) : std::ofstream(name) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  if (!*this) {
    std::cerr << "M_ERROR cannot open file " << name << endl;
    exit(1);
  }
}

M_ofstream::M_ofstream(const char *name, int i) : ofstream(NumberName(name, namebuffer, i)) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  if (!*this) Exit(MError + name)
}

M_ofstream::M_ofstream(const char *name, int i, const char *ext) :
    ofstream(NumberName(name, namebuffer, i, ext)) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  if (!*this) Exit(MError + name)
}

bool FileExists(const char *name) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  std::ifstream file(name);
  if (!file) return false;
  return true;
}

const char *CheckMode(const char *name, const char *mode) {
  if (strcmp(mode, "rename") == 0) {
    if (FileExists(name)) Rename(name);
  } else {
    if (FileExists(name)) Rename(name, mode);
  }
  return name;
}

M_ofstream::M_ofstream(const char *name, const char *mode) : ofstream(CheckMode(name, mode)) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  if (!*this) Exit(MError + name)
}

void M_ofstream::open(const char *name, const char *mode) {
  this->ofstream::open(CheckMode(name, mode));
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  if (!*this) Exit(MError + name)
}

void M_ofstream::open(const char *name) {
  this->ofstream::open(name);
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  if (!*this) Exit(MError + name)
}

void M_ofstream::open_dx(const char *name) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  string Name(name);
  Name += ".dx";
  this->ofstream::open(Name.c_str());
  if (!*this) Exit(MError + name)
}

void M_ofstream::open_gmv(const char *name) {
  Assert(!ParallelProgrammingModel::IsInitialized() || PPM->Master(0));
  string Name(name);
  Name += ".gmv";
  this->ofstream::open(Name.c_str());
  if (!*this) Exit(MError + name)
}

void M_ofstream::popen(const char *name) {
  this->ofstream::open(pNumberName(name, namebuffer));
  if (!*this) Exit(MError + name)
}

void M_ofstream::popen(const char *name, int i) {
  this->ofstream::open(pNumberName(name, namebuffer, i));
  if (!*this) Exit(MError + name)
}

void M_ofstream::popen(const char *name, const char *ext) {
  this->ofstream::open(pNumberName(name, namebuffer, ext));
  if (!*this) Exit(MError + name)
}

void M_ofstream::popen(const char *name, int i, const char *ext) {
  this->ofstream::open(pNumberName(name, namebuffer, i, ext));
  if (!*this) Exit(MError + name)
}
