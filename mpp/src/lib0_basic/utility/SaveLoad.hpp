#ifndef _SAVELOAD_H_
#define _SAVELOAD_H_

#include <complex>
#include <cstdio>
#include <rpc/rpc.h>


#ifdef BUILD_IA

class IAInterval;

class IACInterval;

#endif

/**
 * Provides a Saver and Loader for the basic types. The output ordering has to be the same
 * as the input ordering. Note that no information about the source name / variable is stored
 *
 * Loading and saving operators should be implemented in each class.
 * See e.g. Point::load, Point::save and the corresponding inline operators
 */

bool file_exists(const char *name);

class SaveLoad {
  FILE *file = nullptr;
  XDR xdrs;
public:
  SaveLoad() {}

  ~SaveLoad() { close(); }

  void close();

  static bool FileExists(const char *name, int i = -1);

  static bool PFileExists(const char *name);
protected:
  void open(const char *name, const char *mode, int i = -1);

  void openOnProc(int proc, const char *name, const char *mode, int i = -1);

  XDR *xdr() { return &xdrs; }
};

class Saver : public SaveLoad {
public:
  Saver() : SaveLoad() {}

  explicit Saver(const char *name) : SaveLoad() { open(name); }

  Saver(const char *name, int i) : SaveLoad() { open(name, i); }

  void open(const char *name) { SaveLoad::open(name, "w"); }

  void open(const char *name, int i) { SaveLoad::open(name, "w", i); }

  Saver &operator<<(bool);

  Saver &operator<<(short int);

  Saver &operator<<(int);

  Saver &operator<<(size_t);

  Saver &operator<<(double);

  Saver &operator<<(char);

  Saver &operator<<(std::complex<int> z);

  Saver &operator<<(std::complex<double>);

#ifdef BUILD_IA

  Saver &operator<<(IAInterval);

  Saver &operator<<(IACInterval);

#endif
};

class PSaver : public Saver {
public:
  explicit PSaver(const char *name);
};

// For use on master only
class MSaver : public Saver {
public:
  explicit MSaver(const char *name);
};

class Loader : public SaveLoad {
public:
  /// Checks if files exist on all procs
  static bool files_exist(const char *name);

  Loader() : SaveLoad() {}

  explicit Loader(const char *name) : SaveLoad() { open(name); }

  Loader(const char *name, int i) : SaveLoad() { open(name, i); }

  void open(const char *name) { SaveLoad::open(name, "r"); }

  void open(const char *name, int i) { SaveLoad::open(name, "r", i); }

  Loader &operator>>(bool &);

  Loader &operator>>(short int &);

  Loader &operator>>(int &);

  Loader &operator>>(size_t &);

  Loader &operator>>(double &);

  Loader &operator>>(char &);

  Loader &operator>>(std::complex<double> &);

#ifdef BUILD_IA

  Loader &operator>>(IAInterval &);

  Loader &operator>>(IACInterval &);

#endif
};

class PLoader : public Loader {
public:
  explicit PLoader(const char *name);
};

// For use on master only
class MLoader : public Loader {
public:
  explicit MLoader(const char *name);
};

#endif // of #ifndef _SAVELOAD_H_
