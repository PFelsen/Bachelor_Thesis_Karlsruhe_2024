#ifndef M_IOFILES_HPP
#define M_IOFILES_HPP

#include <fstream>
#include <list>
#include <map>
#include <unordered_map>
#include <vector>

class M_ifstream : public std::ifstream {
public:
  M_ifstream(const char *name, bool test = true);
};

class M_ofstream : public std::ofstream {
public:
  M_ofstream() {}

  M_ofstream(const char *name);

  M_ofstream(const char *, int);

  M_ofstream(const char *, int, const char *);

  M_ofstream(const char *, const char *);

  void open(const char *, const char *);

  void open(const char *);

  void open_dx(const char *);

  void open_gmv(const char *);

  void popen(const char *);

  void popen(const char *, int);

  void popen(const char *, const char *);

  void popen(const char *, int, const char *);
};

bool FileExists(const char *);

inline bool FileExists(const std::string &s) { return FileExists(s.c_str()); };

// inline char *Number_Name(const char *name, int i) {
//     static char NameBuffer[256];
//     return NumberName(name, NameBuffer, i);
// }

#endif // M_IOFILES_HPP
