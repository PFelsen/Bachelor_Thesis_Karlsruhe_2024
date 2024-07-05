#ifndef _CTOOLS_H_
#define _CTOOLS_H_

#include <string>
void Rename(const char *fname, std::string format = "");

char *NumberName(const char *path, const char *name, const char *ext, char *namebuffer, int i);

char *NumberName(const char *name, char *namebuffer, int i);

char *NumberOldName(const char *name, char *namebuffer, int i);

char *NumberName(const char *name, char *namebuffer, int i, const char *ext);

char *pNumberName(const char *name, char *namebuffer);

char *pNumberName(const char *name, char *namebuffer, int i);

char *pNumberName(const char *name, char *namebuffer, int i, const char *ext);

char *pNumberName(const char *name, char *namebuffer, const char *ext);

char *pNumberOldName(const char *name, char *namebuffer);

char *pNumberOldName(const char *name, char *namebuffer, int i);

#endif // of #ifndef _CTOOLS_H_
