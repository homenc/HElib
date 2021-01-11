
class EncodedPtxt;
// container holding either a polynomial
// encoding of either a BGV or CKKS constant
// Interface: not really needed for normal users

class PtxtArray;
// container holding either a either a BGV or CKKS constant
// Interface: see below

class View
{
  // NOTE: View is really a typedef for EncryptedArray
  // but I suggest we start moving away from that

public:
  //========== ENCODING =========
  // BGV-only
  void encode(EncodedPtxt& eptxt, const std::vector<NTL::ZZX>& array) const;
  void encode(EncodedPtxt& eptxt, const std::vector<long>& array) const;

  // CKKS-only encoding functions
  // mag: defaults to Norm(array).
  // prec: defaults to r=getAlMod().getR(), which
  // is usually the same as context.getDefaultPrecision().

  // mag should be an upper bound on Norm(array).
  // If an encoding will be encrypted, the user may wish
  // to hide Norm(array) by setting mag to some data-independent
  // upper bound. A warning is issued if Norm(array) > mag.

  // The encoding will normally have an accuracy of 2^{-prec}, meaning that
  // Norm(array - decode(encode(array))) <= 2^{-prec}.
  // Note that prec may be positive, negative, or zero.
  // The exact logic is a bit heuristic, and a warning is
  // issued if the the accuracy exceeds 2^{-prec}.

  // NOTE: Norm above is the infinity (i.e., max) norm.

  void encode(EncodedPtxt& eptxt,
              const std::vector<cx_double>& array,
              double mag = -1,
              OptLong prec = OptLong()) const;
  void encode(EncodedPtxt& eptxt,
              const std::vector<double>& array,
              double mag = -1,
              OptLong prec = OptLong()) const;

  // BGV and CKKS
  // The following encoding functions are provided for both
  // BGV and CKKS to allow for the creation of "masks".
  // For CKKS, mag and prec are set to their default values.
  void encode(EncodedPtxt& eptxt, const std::vector<bool>& array) const;
  void encodeUnitSelector(EncodedPtxt& eptxt, long i) const;

  //========= ENCRYPTING ==========
  // BGV only
  void encrypt(Ctxt& ctxt, const std::vector<NTL::ZZX>& array) const;
  void encrypt(Ctxt& ctxt, const std::vector<long>& array) const;

  // CKKS only
  // NOTE: mag must be set to non-default value
  void encrypt(Ctxt& ctxt,
               const std::vector<cx_double>& array,
               double mag = -1,
               OptLong prec = OptLong()) const;
  void encrypt(Ctxt& ctxt,
               const std::vector<double>& array,
               double mag = -1,
               OptLong prec = OptLong()) const;

  //========= DECRYPTING ==========
  // BGV-only
  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               std::vector<NTL::ZZX>& ptxt) const;
  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               std::vector<long>& ptxt) const;

  // CKKS-only
  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               std::vector<cx_double>& ptxt,
               OptLong prec = OptLong()) const;
  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               std::vector<double>& ptxt,
               OptLong prec = OptLong()) const;
};

// It is recommended to use the PtxtArray class to do encoding,
// encrypting, and decrypting.
// The interface is a bit more uniform and convenient.

class PtxtArray
{
public:
  // constructors
  explicit PtxtArray(const View& view);
  explicit PtxtArray(const Context& context); // use default view from context
  // NOTE: the View object associated with a PtxtArray
  // object never changes

  PtxtArray(const PtxtArray& other); // copy

  // templates that allow construction via convert:
  // T can be any type supported by convert(PtxtArray,T)
  template <class T>
  PtxtArray(const View& view, const T& t);
  template <class T>
  PtxtArray(const Context& context, const T& t);

  // assignment
  PtxtArray& operator=(const PtxtArray& other); // copy

  // template that allow assignment via load:
  // T can be any type supported by PtxtArray::load(T)
  template <class T>
  PtxtArray& operator=(const T& t);

  void randomReal();
  // CKKS only: random number in [0,1] in each slot

  void randomComplex();
  // CKKS only: random number in complex unit sphere in each slot

  void random();
  // BGV: random ring element in each slot
  // CKKS:  random number in [0,1] in each slot

  // access methods
  const View& getView() const; // preferred name
  const View& getEA() const;   // legacy name

  // The following encode, encrypt, and decrypt operations
  // work for both BGV and CKKS

  void encode(EncodedPtxt& eptxt,
              double mag = -1,
              OptLong prec = OptLong()) const;
  // eptxt = encoding of *this
  // NOTE: for BGV, mag,prec are ignored

  void encrypt(Ctxt& ctxt, double mag = -1, OptLong prec = OptLong()) const;
  // ctxt = encryption of *this
  // NOTES: (1) for BGV, mag,prec are ignored;
  // (2) for CKKS, default mag is set to 2^(ceil(log2(max(Norm(x),1)))),
  // where x is the underlying vector

  void decrypt(const Ctxt& ctxt, const SecKey& sKey, OptLong prec = OptLong());
  // *this = decryption of ctxt under sKey
  // prec is ignored for BGV

  // The following are for CKKS only
  void decryptReal(const Ctxt& ctxt,
                   const SecKey& sKey,
                   OptLong prec = OptLong());
  void decryptComplex(const Ctxt& ctxt,
                      const SecKey& sKey,
                      OptLong prec = OptLong());
  void rawDecrypt(const Ctxt& ctxt, const SecKey& sKey);
  void rawDecryptReal(const Ctxt& ctxt, const SecKey& sKey);
  void rawDecryptComplex(const Ctxt& ctxt, const SecKey& sKey);

  //===== store =====
  // conversion from PtxtArray to std::vector's

  void store(std::vector<NTL::ZZX>& vec) const;  // BGV only
  void store(std::vector<long>& vec) const;      // BGV & CKKS
  void store(std::vector<double>& vec) const;    // CKKS only
  void store(std::vector<cx_double>& vec) const; // CKKS only

  // NOTE: For CKKS, conversion to vector<cx_double> projects
  // the real part, and conversion to vector<long> projects and
  // rounds the real part.

  //===== load =====
  // conversion to PtxtArray from various types
  //(this is a much more permissive set of conversions)

  // vectors
  void load(const std::vector<NTL::ZZX>& vec);  // BGV only
  void load(const std::vector<int>& vec);       // BGV & CKKS
  void load(const std::vector<long>& vec);      // BGV & CKKS
  void load(const std::vector<double>& vec);    // CKKS only
  void load(const std::vector<cx_double>& vec); // CKKS only

  // scalars: puts the same value in each slot
  void load(const NTL::ZZX& scalar); // BGV only
  void load(int scalar);             // BGV & CKKS
  void load(long scalar);            // BGV & CKKS
  void load(double scalar);          // CKKS only
  void load(cx_double scalar);       // CKKS only

  // NTL vectors of NTL ring types
  void load(const NTL::Vec<NTL::GF2>& vec);   // BGV only
  void load(const NTL::Vec<NTL::GF2X>& vec);  // BGV only
  void load(const NTL::Vec<NTL::zz_p>& vec);  // BGV only
  void load(const NTL::Vec<NTL::zz_pX>& vec); // BGV only

  // NTL scalar ring types: puts the same value in each slot
  void load(NTL::GF2 scalar);          // BGV only
  void load(const NTL::GF2X& scalar);  // BGV only
  void load(NTL::zz_p scalar);         // BGV only
  void load(const NTL::zz_pX& scalar); // BGV only

  // NOTE: for conversion to PtxtArray from vectors, the input
  // vector will effectively be truncated or zero-padded to
  // match the number of slots in a PtxtArray
};

std::ostream& operator<<(std::ostream& s, const PtxtArray& a);

bool operator==(const PtxtArray& a, const PtxtArray& b);
bool operator!=(const PtxtArray& a, const PtxtArray& b);
// The above are exact comparisons.
// They are not very useful for CKKS.

// norm and distance functions
double Norm(const PtxtArray& a);
double Distance(const PtxtArray& a, const PtxtArray& b);
// NOTES:
// (1) for BGV, the underlying norm is the trivial norm.
// (2) for CKKS, the underlying norm is the l-infty norm.

// These functions, together with the functions defined in NumbTh.h,
// allow one to write
a ==
    Approx(b)
    // or
    a
    !=
    Approx(b)

    // The Approx function takes two optional parameters:
    double tolerance; // default is 0.01
double floor;         // default is 1.0

// The expression
a == Approx(b, tolerance, floor)
    // is true iff Distance(a,b) <= tolerance*max(Norm(b),floor).
    // The idea is that it checks if the relative error is
    // at most tolerance, unless Norm(b) itself is too small
    // (as determined by floor).
    // For BGV, with the default parameters, a == Approx(b)
    // is equivalent to an exact equality test.

    // arithmetic:
    PtxtArray&
    operator+=(PtxtArray& a, const PtxtArray& b);
PtxtArray& operator-=(PtxtArray& a, const PtxtArray& b);
PtxtArray& operator*=(PtxtArray& a, const PtxtArray& b);

// also for b of of any type T supported by PtxtArray::load(T):
// NOTE: SFINAE is used to truly limit T
template <class T>
PtxtArray& operator+=(PtxtArray& a, const T& b);
template <class T>
PtxtArray& operator-=(PtxtArray& a, const T& b);
template <class T>
PtxtArray& operator*=(PtxtArray& a, const T& b);

void negate(PtxtArray& a);
void power(PtxtArray& a, long e);

// data movement
void rotate(PtxtArray& a, long k);
void shift(PtxtArray& a, long k);

voud rotate1D(PtxtArray& a, long i, long k);
voud shift1D(PtxtArray& a, long i, long k);

void applyPerm(PtxtArray& a, const NTL::Vec<long>& pi);

// For CKKS, these use complex conjugation, rather than Frobenius
void frobeniusAutomorph(PtxtArray& a, long j);
void frobeniusAutomorph(PtxtArray& a, const NTL::Vec<long>& vec);
void conjugate(PtxtArray& a); // synonym for frobeniusAutomorph(ctxt, 1);

void extractRealPart(PtxtArray& a); // CKKS only
void extractImPart(PtxtArray& a);   // CKKS only

void totalSums(PtxtArray& a);
void runningSums(PtxtArray& a);
