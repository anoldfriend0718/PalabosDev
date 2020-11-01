#pragma once

#include "parameters/LBMModelParser2D.h"
namespace plb {
namespace util {

template <typename T_, template <typename U> class Descriptor>
std::string getDynamicsName(int id) {
  plb::meta::DynamicsRegistration<T_, Descriptor> &dynamicsRegistration =
      plb::meta::dynamicsRegistration<T_, Descriptor>();
  std::string dynamicsName = dynamicsRegistration.getName(id);
  return dynamicsName;
}

template <typename T_, template <typename U> class Descriptor>
std::string getDynamicsName(LBMModelParser2D<T_, Descriptor> &model) {
  int id = model.getDynamics()->getId();
  std::string dynamicsName = getDynamicsName<T_, Descriptor>(id);
  return dynamicsName;
}

std::string getFileName(const std::string path,
                        const std::string seperator = "/") {
  std::string::size_type iPos = path.find_last_of(seperator) + 1;
  std::string filename = path.substr(iPos, path.length() - iPos);
  std::string name = filename.substr(0, filename.rfind("."));
  return name;
}

} // namespace util
} // namespace plb