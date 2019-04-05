//Header guard
#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <exception>
#include <stdexcept>
#include <sstream>

namespace helib {

class Exception
{
  public:
    virtual ~Exception() = default;
    virtual const char* what() const noexcept = 0;
  protected:
    Exception() = default;
};
  
class LogicError : public std::logic_error, public ::helib::Exception
{
public:
  explicit LogicError(const std::string& what_arg) : std::logic_error(what_arg) {};
  explicit LogicError(const char* what_arg) : std::logic_error(what_arg) {};
  virtual ~LogicError(){};
  virtual const char* what() const noexcept override {return std::logic_error::what();};
};

// Done up to here

class OutOfRangeError : public std::out_of_range, public ::helib::Exception
{
public:
  explicit OutOfRangeError(const std::string& what_arg) : std::out_of_range(what_arg) {};
  explicit OutOfRangeError(const char* what_arg) : std::out_of_range(what_arg) {};
  virtual ~OutOfRangeError(){};
  virtual const char* what() const noexcept override {return std::out_of_range::what();};
};

class RuntimeError : public std::runtime_error, public ::helib::Exception
{
public:
  explicit RuntimeError(const std::string& what_arg) : std::runtime_error(what_arg) {};
  explicit RuntimeError(const char* what_arg) : std::runtime_error(what_arg) {};
  virtual ~RuntimeError(){};
  virtual const char* what() const noexcept override {return std::runtime_error::what();};
};
  
class InvalidArgument : public std::invalid_argument, public ::helib::Exception
{
public:
  explicit InvalidArgument(const std::string& what_arg) : std::invalid_argument(what_arg) {};
  explicit InvalidArgument(const char* what_arg) : std::invalid_argument(what_arg) {};
  virtual ~InvalidArgument(){};
  virtual const char* what() const noexcept override {return std::invalid_argument::what();};
};
  
  
} // End of namespace
#endif // End of header guard
