#include <sys/stat.h>

#include "Parallel.hpp"

void Rename(const char *fname, std::string format) {
  std::string new_fname(fname);
  if (format != "") {
    format = "." + format;
    if (new_fname.substr(new_fname.length() - format.size(), format.size()) == format) {
      new_fname.erase(new_fname.length() - format.size(), format.size());
    }
  }
  new_fname.append("_");

  struct stat fstat;
  if (stat(fname, &fstat) < 0) exit(1);
  time_t T = fstat.st_mtime;
  char timestamp[128];

  strftime(timestamp, 64, "%y-%m-%d_%H-%M-%S", localtime(&T));
  new_fname.append(timestamp);
  new_fname.append(format);
  if (rename(fname, new_fname.c_str())) exit(1);
}

char *NumberName(const char *path, const char *name, const char *ext, char *namebuffer, int i) {
  if (i < 0) sprintf(namebuffer, "%s/%s.%s", path, name, ext);
  else sprintf(namebuffer, "%s/%s.%04d.%s", path, name, i, ext);
  return namebuffer;
}

char *NumberName(const char *name, char *namebuffer, int i) {
  if (i < 0) sprintf(namebuffer, "%s", name);
  else sprintf(namebuffer, "%s.%04d", name, i);
  return namebuffer;
}

char *NumberOldName(const char *name, char *namebuffer, int i) {
  if (i < 0) sprintf(namebuffer, "%s.old", name);
  else sprintf(namebuffer, "%s.old.%04d", name, i);
  return namebuffer;
}

char *NumberName(const char *name, char *namebuffer, int i, const char *ext) {
  if (i < 0) sprintf(namebuffer, "%s.%s", name, ext);
  else sprintf(namebuffer, "%s.%04d.%s", name, i, ext);
  return namebuffer;
}

char *pNumberName(const char *name, char *namebuffer) {
  sprintf(namebuffer, "%s.p%04d", name, PPM->Proc(0));
  return namebuffer;
}

char *pNumberName(const char *name, char *namebuffer, int i) {
  if (i < 0) sprintf(namebuffer, "%s.p%04d", name, PPM->Proc(0));
  else sprintf(namebuffer, "%s.p%04d.%04d", name, PPM->Proc(0), i);
  return namebuffer;
}

char *pNumberName(const char *name, char *namebuffer, int i, const char *ext) {
  sprintf(namebuffer, "%s.p%04d.%04d.%s", name, PPM->Proc(0), i, ext);
  return namebuffer;
}

char *pNumberName(const char *name, char *namebuffer, const char *ext) {
  sprintf(namebuffer, "%s.p%04d.%s", name, PPM->Proc(0), ext);
  return namebuffer;
}

char *pNumberOldName(const char *name, char *namebuffer) {
  sprintf(namebuffer, "%s.old.p%04d", name, PPM->Proc(0));
  return namebuffer;
}

char *pNumberOldName(const char *name, char *namebuffer, int i) {
  if (i < 0) sprintf(namebuffer, "%s.old.p%04d", name, PPM->Proc(0));
  sprintf(namebuffer, "%s.old.p%04d.%04d", name, PPM->Proc(0), i);
  return namebuffer;
}