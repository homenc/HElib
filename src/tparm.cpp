namespace std {} using namespace std;
namespace NTL {} using namespace NTL;

#include "NumbTh.h"

#if 0

/* doArgProcessing: converts c-string s to value T,
 * returns upon success.  By default, we parse using
 * the istream input operator, except when T = string
 * and just convert without any parsing.
 */

template<class T>
bool doArgProcessing(T *value, const char *s)
{
  string ss(s);
  stringstream sss(ss);
  return sss >> *value;
}

bool doArgProcessing(string *value, const char *s)
{
  *value = string(s);
  return true;
}

/* ArgProcessor: virtual base class */

class ArgProcessor {
public:
virtual bool process(const char *s) = 0;
};

/* ArgProcessorDerived: templated subclasses */

template<class T>
class ArgProcessorDerived : public ArgProcessor   {
public:
  T *value;

  virtual bool process(const char *s)
  {
    return doArgProcessing(value, s);
  }

  ArgProcessorDerived(T* _value) : value(_value) {}
};

class ArgMapping {
public:
  unordered_map< string, shared_ptr<ArgProcessor> > map;
  stringstream doc;

  template<class T>
  void arg(const char *name, T& value) 
  { 
    shared_ptr<ArgProcessor> ap = 
      shared_ptr<ArgProcessor>(new ArgProcessorDerived<T>(&value));

    map[name] = ap;
  }

  template<class T>
  void arg(const char *name, T& value, const char *doc1) 
  {
    arg(name, value);
    doc << "\t" << name << " \t" << doc1 << " \t[" << value << "]" << "\n";
  }


  bool parse(int argc, const char *argv[])
  {
    for (long i = 1; i < argc; i++) {
      const char *x = argv[i];
      long j = 0;
      while (x[j] != '=' && x[j] != '\0') j++; 
      if (x[j] == '\0') return false;
      string name(x, j);
      const char *s = x+j+1;

      shared_ptr<ArgProcessor> ap = map[name];
      if (!ap) return false;
      if (!ap->process(s)) return false;
    }

    return true;
  }

  string documentation() 
  {
    return doc.str();
  }
};

#endif

int main(int argc, const char *argv[])
{
  ArgMapping amap;

  long p = 2;
  amap.arg("p", p, "the number p");

  double v = 1.2;
  amap.arg("v", v, "the vector v");

  if (amap.parse(argc, argv))
    cout << "OK\n";
  else
    cout << "BAD\n";

  cout << p << "\n";
  cout << v << "\n";

  cout << amap.documentation();

}
