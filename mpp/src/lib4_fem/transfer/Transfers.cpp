#include "Transfers.hpp"
#include "LagrangeTransfer.hpp"
#include "MatrixTransfer.hpp"
#include "TransferWrapper.hpp"

std::unique_ptr<ITransfer> GetTransfer(const Vector &coarse, const Vector &fine, std::string name) {
  const auto &coarseType = typeid(coarse.GetDisc());
  const auto &fineType = typeid(fine.GetDisc());

  if (coarseType != fineType) { THROW("Discretizations of Vectors don't match") }

  Config::Get("Transfer", name);
  if (contains(name, "Matrix")) {
    return std::make_unique<ITransfer>(new MatrixTransfer(coarse, fine));
  }

  if (contains(name, "Legacy")) {
    return std::make_unique<ITransfer>(new TransferWrapper(name, coarse, fine));
  }

  if (fineType == typeid(LagrangeDiscretization)) {
    return std::make_unique<ITransfer>(new LagrangeTransfer(coarse, fine));
  }

  return std::make_unique<ITransfer>(new TransferWrapper(name, coarse, fine));
}
