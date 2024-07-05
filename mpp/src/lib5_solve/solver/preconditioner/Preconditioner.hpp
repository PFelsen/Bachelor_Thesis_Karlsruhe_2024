#ifndef _PRECONDITIONER_H_
#define _PRECONDITIONER_H_

#include "Assemble.hpp"
#include "PreconditionerBase.hpp"
#include "Sparse.hpp"
#include "TimeDate.hpp"

std::string PCName(const string &prefix = "");

Preconditioner *GetPC();

Preconditioner *GetPC(const string &);

Preconditioner *GetPCByPrefix(const string &);

#endif // of #ifndef _PRECONDITIONER_H_
