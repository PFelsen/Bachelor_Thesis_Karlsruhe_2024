#include "ProblemBase.hpp"


std::string to_string(COMPONENT comp) {
  switch (comp) {
    case COMPONENT::V_X: return "V_X";
    case COMPONENT::V_Y: return "V_Y";
    case COMPONENT::V_Z: return "V_Z";
    case COMPONENT::P0: return "P0";
    case COMPONENT::P1: return "P1";
    case COMPONENT::P2: return "P2";
    case COMPONENT::P3: return "P3";
    case COMPONENT::P4: return "P4";
    case COMPONENT::P5: return "P5";
    default:
      Exit("Implement std::to_string(Component)")
  }
}

bool isVelocity(COMPONENT comp) {
  switch (comp) {
    case COMPONENT::V_X:
    case COMPONENT::V_Y:
    case COMPONENT::V_Z: return true;
    case COMPONENT::P0:
    case COMPONENT::P1:
    case COMPONENT::P2:
    case COMPONENT::P3:
    case COMPONENT::P4:
    case COMPONENT::P5: return false;
    default:
    Exit("Implement isVelocity(Component)")
  }
}

int toVelocityIndex(COMPONENT comp) {
  switch (comp) {
    case COMPONENT::V_X: return 0;
    case COMPONENT::V_Y: return 1;
    case COMPONENT::V_Z: return 2;
    default:
    Exit("Errpr called std::toVelocityIndex(Component)")
  }
}

COMPONENT GetDampingComponent(int i) {
  switch (i) {
    case 0:return COMPONENT::P1;
    case 1:return COMPONENT::P2;
    case 2:return COMPONENT::P3;
    case 3:return COMPONENT::P4;
    case 4:return COMPONENT::P5;
    default: ERROR("Damping not implemented for i=" + std::to_string(i))
  }
}

double GetDampingIndex(COMPONENT comp) {
  switch (comp) {
    case COMPONENT::P1:return 0;
    case COMPONENT::P2:return 1;
    case COMPONENT::P3:return 2;
    case COMPONENT::P4:return 3;
    case COMPONENT::P5:return 4;
    default: ERROR("Damping not implemented for comp")
  }
}