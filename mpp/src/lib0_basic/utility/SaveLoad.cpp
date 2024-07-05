#include <complex>
#include <fstream>

#include "Parallel.hpp"
#include "SaveLoad.hpp"
#include "ctools.hpp"


#ifdef BUILD_IA

#include "IACInterval.hpp"

#endif

bool file_exists(const char *name) {
  std::ifstream file(name);
  if (!file) return false;
  return true;
}

void SaveLoad::open(const char *name, const char *mode, int i) {
  char namebuffer[128];
  pNumberName(name, namebuffer, i);
  if (mode[0] == 'w') {
    if (file_exists(namebuffer)) {
      char old[128];
      if (rename(namebuffer, pNumberOldName(name, old, i))) { THROW("Error: Cannot copy files") }
    }
  }
  file = fopen(pNumberName(name, namebuffer, i), mode);
  if (PPM->Or(!file, 0)) THROW("Error: Cannot open files")

  if (mode[0] == 'w') {
    xdrstdio_create(&xdrs, file, XDR_ENCODE);
    int procs = PPM->Size(0);
    xdr_int(xdr(), &procs);
  }
  if (mode[0] == 'r') {
    xdrstdio_create(&xdrs, file, XDR_DECODE);
    int procs;
    xdr_int(xdr(), &procs);
    if (procs != int(PPM->Size(0))) { WarningOnProc("Error: Number of procs does not fit") }
  }
}

void SaveLoad::openOnProc(int proc, const char *name, const char *mode, int i) {
  if (PPM->Proc(0) != proc) return;

  char namebuffer[128];
  NumberName(name, namebuffer, i);
  if (mode[0] == 'w') {
    if (file_exists(namebuffer)) {
      char old[128];
      if (rename(namebuffer, NumberOldName(name, old, i))) { THROW("Error: Cannot copy files") }
    }
  }
  file = fopen(NumberName(name, namebuffer, i), mode);

  if (mode[0] == 'w') {
    xdrstdio_create(&xdrs, file, XDR_ENCODE);
    xdr_int(xdr(), &proc);
  }
  if (mode[0] == 'r') {
    xdrstdio_create(&xdrs, file, XDR_DECODE);
    int p;
    xdr_int(xdr(), &p);
    if (p != proc) { WarningOnProc("Error: Number of proc does not fit") }
  }
}

void SaveLoad::close() {
  if (!file) return;
  fclose(file);
  file = nullptr;
}

bool SaveLoad::FileExists(const char *name, int i) {
  char namebuffer[128];
  pNumberName(name, namebuffer, i);
  return PPM->And(file_exists(namebuffer), 0);
}

bool SaveLoad::PFileExists(const char *name) { return FileExists(name, PPM->Size(0)); }

Saver &Saver::operator<<(bool i) {
  bool_t tmp = i;
  xdr_bool(xdr(), &tmp);
  return *this;
}

Saver &Saver::operator<<(short int i) {
  xdr_short(xdr(), &i);
  return *this;
}

Saver &Saver::operator<<(int i) {
  xdr_int(xdr(), &i);
  return *this;
}

Saver &Saver::operator<<(size_t i) {
  int j = int(i); // CONVERSION!
  xdr_int(xdr(), &j);
  return *this;
}

Saver &Saver::operator<<(double a) {
  xdr_double(xdr(), &a);
  return *this;
}

Saver &Saver::operator<<(char a) {
  xdr_char(xdr(), &a);
  return *this;
}

Saver &Saver::operator<<(std::complex<int> z) {
  int re = std::real(z);
  int im = std::imag(z);
  xdr_int(xdr(), &re);
  xdr_int(xdr(), &im);
  return *this;
}

Saver &Saver::operator<<(std::complex<double> z) {
  double re = std::real(z);
  double im = std::imag(z);
  xdr_double(xdr(), &re);
  xdr_double(xdr(), &im);
  return *this;
}

#ifdef BUILD_IA

Saver &Saver::operator<<(IAInterval a) { return *this << inf(a) << sup(a); }

Saver &Saver::operator<<(IACInterval a) { return *this << real(a) << imag(a); }

#endif

PSaver::PSaver(const char *name) : Saver(name, PPM->Size(0)) {}

MSaver::MSaver(const char *name) { SaveLoad::openOnProc(0, name, "w"); }

bool Loader::files_exist(const char *name) {
  char namebuffer[128];
  pNumberName(name, namebuffer, PPM->Size(0));
  return PPM->And(file_exists(namebuffer), 0);
}

Loader &Loader::operator>>(bool &i) {
  bool_t tmp;
  xdr_bool(xdr(), &tmp);
  i = tmp;
  return *this;
}

Loader &Loader::operator>>(short int &i) {
  xdr_short(xdr(), &i);
  return *this;
}

Loader &Loader::operator>>(int &i) {
  xdr_int(xdr(), &i);
  return *this;
}

Loader &Loader::operator>>(size_t &i) {
  int j;
  xdr_int(xdr(), &j);
  i = size_t(j);
  return *this;
}

Loader &Loader::operator>>(double &a) {
  xdr_double(xdr(), &a);
  return *this;
}

Loader &Loader::operator>>(char &a) {
  xdr_char(xdr(), &a);
  return *this;
}

Loader &Loader::operator>>(std::complex<double> &z) {
  double re, im;
  xdr_double(xdr(), &re);
  xdr_double(xdr(), &im);
  z = std::complex<double>(re, im);
  return *this;
}

#ifdef BUILD_IA

Loader &Loader::operator>>(IAInterval &a) {
  double inf, sup;
  *this >> inf >> sup;
  a = IAInterval(inf, sup);
  return *this;
}

Loader &Loader::operator>>(IACInterval &a) {
  IAInterval re, im;
  *this >> re >> im;
  a = IACInterval(re, im);
  return *this;
}

#endif

PLoader::PLoader(const char *name) : Loader(name, PPM->Size(0)) {}

MLoader::MLoader(const char *name) { SaveLoad::openOnProc(0, name, "r"); }