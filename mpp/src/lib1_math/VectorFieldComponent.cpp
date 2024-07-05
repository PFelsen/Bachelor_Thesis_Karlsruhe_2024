#include "VectorFieldComponent.hpp"

#ifdef BUILD_IA

VelocityComponentT<double> mid(const VelocityComponentT<IAInterval> &IAx) {
  return VelocityComponentT<double>(IAx.Component(), mid(IAx.Value()));
}

VelocityComponentT<double> sup(const VelocityComponentT<IAInterval> &IAx) {
  return VelocityComponentT<double>(IAx.Component(), sup(IAx.Value()));
}

VelocityComponentT<double> inf(const VelocityComponentT<IAInterval> &IAx) {
  return VelocityComponentT<double>(IAx.Component(), inf(IAx.Value()));
}


#endif // BUILD_IA
