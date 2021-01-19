#include <helib/JsonWrapper.h>
#include <helib/assertions.h>
#include <string>

#include "io.h"

#include <json.hpp>
using json = ::nlohmann::json;

namespace helib {

std::string JsonWrapper::pretty(long indent) const
{
  try {
    auto j = std::any_cast<json>(this->json_obj);
    return j.dump(indent);
  } catch (const std::bad_any_cast& e) {
    throw LogicError(std::string("Bad cast to a JSON object: \n\t") + e.what());
  }
}

std::string JsonWrapper::toString() const
{
  try {
    auto j = std::any_cast<json>(this->json_obj);
    return j.dump();
  } catch (const std::bad_any_cast& e) {
    throw LogicError(std::string("Bad cast to a JSON object: \n\t") + e.what());
  }
}

JsonWrapper JsonWrapper::at(const std::string& key) const
{
  try {
    auto j = std::any_cast<json>(this->json_obj);
    return wrap(j.at(key));
  } catch (const std::bad_any_cast& e) {
    throw LogicError(std::string("Bad cast to a JSON object: \n\t") + e.what());
  }
}

std::ostream& operator<<(std::ostream& str, const JsonWrapper& wrapper)
{
  return str << wrapper.toString();
}

} // namespace helib
