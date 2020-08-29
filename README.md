HElib
=====

HElib is an open-source ([Apache License v2.0][5]) software library that 
implements [homomorphic encryption][6] (HE). Currently available schemes 
are the implementations of the [Brakerski-Gentry-Vaikuntanathan][1] (BGV) 
scheme with [bootstrapping][8] and the Approximate Number scheme of 
[Cheon-Kim-Kim-Song][9] (CKKS), along with many optimizations to make 
homomorphic evaluation run faster, focusing mostly on effective use of 
the [Smart-Vercauteren][2] ciphertext packing techniques and
the [Gentry-Halevi-Smart][3] optimizations. See [this report][7] for a
description of a few of the algorithms using in this library. 

Since mid-2018 HElib has been under extensive refactoring for *Reliability*, 
*Robustness & Serviceability*, *Performance*, and most importantly *Usability* 
for researchers and developers working on HE and its uses.

HElib supports an *"assembly language for HE"*, providing low-level routines
(set, add, multiply, shift, etc.), sophisticated automatic noise management,
improved BGV bootstrapping, multi-threading, and also support for Ptxt (plaintext) 
objects which mimics the functionality of Ctxt (ciphertext) objects. 
See [changes.md](changes.md) for more details.

Full installation instructions and a list of the required dependencies can be found 
in [INSTALL.md](INSTALL.md).

If you are interested in contributing to HElib, please read our 
[Contributing Guidelines](CONTRIBUTING.md).

HElib is written in C++14 and uses the [NTL mathematical library][4].  
HElib is distributed under the terms of the [Apache License v2.0][5].  


  [1]: http://eprint.iacr.org/2011/277       "BGV12"
  [2]: http://eprint.iacr.org/2011/133       "SV11"
  [3]: http://eprint.iacr.org/2012/099       "GHS12"
  [4]: http://www.shoup.net/ntl/             "NTL"
  [5]: http://www.apache.org/licenses/LICENSE-2.0  "Apache-v2.0"
  [6]: http://en.wikipedia.org/wiki/Homomorphic_encryption "Homomorphic encryption"
  [7]: http://eprint.iacr.org/2014/106       "algorithms"
  [8]: http://eprint.iacr.org/2014/873       "bootstrapping"
  [9]: http://eprint.iacr.org/2016/421       "CKKS16"
  [10]: https://github.com/homenc/HElib      "GitHubPages"
