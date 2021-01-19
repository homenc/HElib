#ifndef JSONWRAPPER_HIDDEN_H
#define JSONWRAPPER_HIDDEN_H
#include <string>
#include <any>

namespace helib {

struct JsonWrapper
{
public:
  explicit JsonWrapper(const std::any& json_repr) : json_obj(json_repr) {}
  explicit operator bool() const { return json_obj.has_value(); }
  const std::any& getJSONobj() const { return json_obj; }
  JsonWrapper at(const std::string& key) const;
  std::string toString() const;

  friend std::ostream& operator<<(std::ostream& str,
                                  const JsonWrapper& wrapper);

  std::string pretty(long indent = 2) const;

private:
  std::any json_obj;
};

} // namespace helib

#endif // JSONWRAPPER_HIDDEN_H
