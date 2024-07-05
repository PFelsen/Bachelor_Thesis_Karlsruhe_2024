#ifndef TRANSFERS_HPP
#define TRANSFERS_HPP

#include <Config.hpp>
#include "ITransfer.hpp"

std::unique_ptr<ITransfer> GetTransfer(const Vector &coarse, const Vector &fine,
                                       std::string name = "");

#endif // TRANSFERS_HPP
