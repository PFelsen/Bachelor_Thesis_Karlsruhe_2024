#include "SpaceTimeTransfer.hpp"


std::unique_ptr<SpaceTimeTransfer> GetSpaceTimeTransfer(const STDiscretization &disc,
                                                        const std::string &name) {
  if (disc.isDgInTime()) {
    if (name == "Projection") {
      return std::make_unique<SpaceTimeTransferProjectionDGDG>(disc);
    } else if (name == "Interpolation") {
      return std::make_unique<SpaceTimeTransferInterpolateDGDG>(disc);
    }
    return std::make_unique<SpaceTimeTransferInterpolateDGDG>(disc);
  }else{
    if (name == "Interpolation"){
      return std::make_unique<SpaceTimeTransferInterpolation>(disc);
    }
  }
  THROW("Transfer with name " + name + " not implemented!");
}