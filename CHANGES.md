HElib 1.1.0, August 2020
===============================
(tagged as v1.1.0)

June-August 2020
------------------
* Utils
  * create-context
  * encrypt/decrypt
  * encode/decode examples
* Improvements to noise management
* Move to CMake >= 3.10.2
* Doxygen improved look and feel
* Improved docs
* New examples and mini-tutorial on private preserving db query
* Configurable PSI/PIR beta
* Some tools in `misc/` for code formatting and algebra generation
* `Warning` messages are now logged (default `helib.log`)
* Bug fixes

HElib 1.0.2, June 2020
===============================
(tagged as v1.0.2)

April-June 2020
------------------
* Source code contribution policy
* HElib now builds with compiler pedantic flags
* Consistent C++ format style adopted
* Improved sampling functions
* Bug fixes

HElib 1.0.1, April 2020
===============================
(tagged as v1.0.1)

January-April 2020
------------------
* Additional BGV examples BGV_(packed and binary)_arithmetic and simple BGV_database_lookup
* Improved powerful basis
* Enhancements and extensions to Binary operations
* Enhancements to Ptxt functions
* Enhancements to Ctxt Functions
* Enhancements and fixes to `ArgMap` functionality
* Compilation warning messages clean-up
* Updates to Docker build
* Using stable Google Test framework to HElib testing suite.
* Additional tests ported to using Google Test framework
* Bug fixes

HElib 1.0.0, January 2020
===============================
(tagged as v1.0.0)

December 2019
-------------
* C++14 Standard (minimum level)
* New `Ptxt` Plaintext class that implements the same functionality of the `Ctxt` ciphertext class.  
* Improved version of the `ArgMap` API for command line arguments.  
* Restructuring of the project directory tree.  
* Removed AES example - improved version on its way.  
* Doxygen documentation.  
* Bug fixes.

HElib 1.0.0 beta 2, November 2019
===============================
(tagged as 1.0.0-beta2-Nov2019)

September-November 2019
-----------------------
* Significant refactoring and cleanup of codebase.  
* New `helib` namespace.  
* New `examples` and `benchmarks` directory trees.  
* Improvements to bootstrapping.  
* Better tests for bootstrapping and binary arithmetic in BGV.   
* Docs and example code for binary arithmetic.  
* Overall code and performance improvements in `NumbTh.cpp`.   
* HElib now avoids *very bad* generators.   
* Bug fixes.

HElib 1.0.0 beta 1, August 2019
===============================
(tagged as 1.0.0-beta1-Aug2019)

August 2019
-----------
* Improved noise management in HElib.  
* Better and more robust bootstrapping algorithm.

July 2019
---------
* Added new bootstrapping and PGFFT tests.

June 2019
---------
* Added implementation of PGFFT to replace Armadillo for FFTs in
the CKKS scheme. See comments in `PGFFT.h` for more information.

May 2019
--------
* CKKS bug fixes and fixed compiler warnings.

April 2019
----------
* Moved the HElib repository on github from shaih/HElib to HomEnc/HElib.  
* Introduced HElib-specific exceptions, replaced C `assert`, NTL
* Error and `std::exception` to use HElib-specific assertions that throw the new exceptions.  
* See comments in `assertions.h` and `exceptions.h` for usage information.

March 2019
----------
* Introduced new test framework (google test), documented in [TESTS.md](TESTS.md).  
* Previous framework will be deprecated. Added a cmake build script for building HElib and dependencies, documented in [INSTALL.md](INSTALL.md).  
* Added an example program, see `example program`.


HElib 1.0.0 beta 0, January 2019
===============================
(tagged as 1.0.0-beta0-Jan2019)
 
This commit includes multiple changes to the library. We hope to add documentation over the next few weeks, then make it into an official version 1.0.0 release. Below is a summary of the main changes since March of 2018.


## 1. Decoupling modulus-switching resolution from size of primes in the chain

With this version, we greatly refined the supported resolution for modulus switching. In previous versions all the modulus in the chain were of the same size, except perhaps one that was exactly half the size of the others. As a result, modulus switching was always confined to a rather coarse resolution. For example when the moduli size was 50 bits and the half-size modulus was 25 bits long, modulus switching had to drop multiples of 25 bits every time.

In the current version the application can specify the target resolution, using the optional 'resolution' argument of the method `FHEcontext::buildModChain` (with `resolution=3` as the default). HElib now keeps moduli of different sizes, and attempts to find a subset of them with sizes adding up to what's needed.

An important implication of this change is that HElib no longer has a notion of a discrete number of "levels" that are available for computation. Instead it has a notion of "capacity" of a ciphertext, which is just the logarithm of the modulus/noise ratio. One can get the current capacity of a ciphertext by calling `ct.capacity()` for the natural logarithm or `ct.bitCapacity()` for the logarithm base 2.

When calling `FHEcontext::buildModChain`, the number-of-levels argument is replaced by a number-of-bits arguments. The task of determining the required number of bits for any given computation is left to the calling application. (We may provide some functions to help make that decisions, but they are not available yet.)

Some convenience methods that we do provide are `ct.naturalSize()` that returns the size (natural logarithm) that this ciphertext should be mod-switched to if we square it, and `ct.naturalPrimeSet()` that returns the corresponding prime-set.


## 2. Preliminary support for the CKKS cryptosystem

Included in this version is a draft implementation of the CKKS approximate-number homomorphic encryption. To access this implementation, one needs to create an `FHEcontext` with "plaintext space modulus" `p=-1`, and using the parameter `r` to specify the requested precision. For example setting
```c++
  FHEcontext context(/*m=*/4096, /*p=*/-1, /*r=*/8);
```
will set up an instance of CKKS over the 2^12 cyclotomic ring, requesting to keep precision of 1/256 (i.e., 8 bits to the right of the binary point). The implementation includes a `EncryptedArrayCx` class that handles the encoding and decoding of these instances, see EncryptedArray.h.

We use the same chassis for CKKS as for BGV, and in particular we support arbitrary cyclotomics (not just powers of two). But not everything is implemented for this case yet. (For instance `EncryptedArray::shift1D` is not implemented for CKKS ciphertexts.)

This implementation is quite preliminary, and for now it only "does the right thing" as long as all the complex data elements throughout the computation are close to one in absolute value. For smaller data values the requested precision will typically not be enough, while for larger values the implementation will spend too much resources trying to keep the precision way too high.


## 3. More rigorous dynamic noise estimates

HElib includes logic to estimate the level of noise in the ciphertext, and many of the library's housekeeping operations depend on this estimate. Crucially, big discrepancies between the actual noise and the estimated one may lead to decryption errors.

Earlier version roughly kept an estimate of the L2-squared noise of the coefficient representation of the ring elements, using various hand-wavy methods (which often times were not even heuristically sound). The current version instead keeps high-probability bound on the L-infinity norm in the canonical embedding, and tries to provide at least a heuristic guarantee of no decryption errors. The result is sometimes more conservative modulus-switching choices.

We added to the `FHEcontext` class a data member `scale` that represents the "number of standard deviations" used in our high-probability heuristic bounds. It defaults to `scale=10`, corresponding to probability roughly 2^{-76} for a normal random variable being more than 10 standard deviations away from its mean. The application can set `context.scale` to other values, it controls the estimated size of the noise terms in freshly sampled elements and other operations (such as modulus switching and key switching). See more documentation in FHEContext.h.

Also with this version we no longer use the "ring constant" data member `cM` in the PAlgebra class for noise estimate. That data member is only used in recryption, to tweak the size computation for the "powerful basis" of elements. See more documentation in recryption.h under the `setAE()` function.


## 4. Other notable changes

### 4a. Keeping an integer factor in BGV
In previous versions we used the decryption invariant

    [<sk,ct>]_q = q*m (mod p),

which is convenient for the modulus switching operation. One of the consequences of this invariant, however, is that whenever we multiply two ciphertexts we have to multiply the result by q^{-1} mod p, which increases the noise (especially for large plaintext spaces p).

In this version we added a `long intFactor` data member to ciphertexts, and we now use the modified invariant

    [<sk,ct>]_q = intFactor*q*m (mod p),

This is just as convenient as before for modulus switching, but it allows us to just modify the `intFactor` (without changing the noise) after multiplication. When we add two ciphertexts we still need to make sure that they both have the same `intFactor`, but it is easy to see that this can be done while increasing the noise of the result by at most a sqrt(p) factor.

### 4b. Wider noise sampling for non-power-of-two ring-LWE
Fixed a security bug, related to the ring-LWE assumption in non-power-of-two cyclotomic rings. Before we always sampled the noise with a constant width in the coefficient representation (sigma=3.2 by default). This is an acceptable choice for power-of-two cyclotomics, but not otherwise. In the new version, for the m'th cyclotomic (with m not a power of two), we sample a *degree-m polynomial* in coefficient representation using Gaussian width sigma*sqrt(m), and then reduce the result modulo Phi_m(X). This yields somewhat larger noise terms.

On the other hand, when sampling noise terms during key-generation, we check the canonical-embedding norm of the result and re-sample if it is too large. (Specifically we set the parameters so that the probability of re-sampling is below 1/2.) This very often yields smaller noise terms for the keys.

### 4c. Implemented the Chen-Han "thin" bootstrapping procedure
Implemented the faster procedure for bootstrapping lightly packed ciphertexts (where the slots contain only integers), from [[Chen-Han, Eurocrypt 2018]](https://ia.cr/2018/067). Integrated this procedure with the improved linear algebra methods of HElib.

### 4d. Bootstrapping parameters handled more rigorously
Corrected some mistakes and tightened the analysis in [[Halevi-Shoup, Eurocrypt 2015]](https://ia.cr/2014/873), specifically Lemmas 5.1 and 5.2 and appendix A. This makes very little difference for bootstrapping with a prime plaintext space p, but a dramatic difference for plaintext space p^r for r>1. For example bootstrapping with plaintext space p^r=2^8 is now almost as efficient as for plaintext space p=2.

## 5. Minor changes
* Modified `Makefile` to look for local settings in the file `local-defs`;
* Introduced binary serialization/deserialization (vs. the previous ascii-based method);
* Removed all 'using' from header files;
* Added innerProduct function on wrapped vectors (see `CtPtrs.h`);
* `DoubleCRT` constructors no longer have default prime-set;
* Added a "destructive" multiply method, calling `ct1.multLowLvl(ct, true)` may mod-switch both `c1` and `c2`;
* Made `FHEPubKey::Encrypt` virtual, overridden by `FHESecKey::Encrypt`;
* Small changes in how the digits are set-up in `FHEcontext` to make their size more uniform;
* Eliminated some code duplication between "thin" and "think" bootstrapping in recryption.cpp;
* Added some more functions/options in the debugging module.
