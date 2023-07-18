HElib
=====

[![Build and Test](https://github.com/homenc/HElib/actions/workflows/github-ci.yml/badge.svg)](https://github.com/homenc/HElib/actions/workflows/github-ci.yml)

HElib is an open-source ([Apache License v2.0][5]) software library that
implements [homomorphic encryption][6] (HE). Currently available schemes are the
implementations of the [Brakerski-Gentry-Vaikuntanathan][1] (BGV) scheme with
[bootstrapping][8] and the Approximate Number scheme of [Cheon-Kim-Kim-Song][9]
(CKKS), along with many optimizations to make homomorphic evaluation run faster,
focusing mostly on effective use of the [Smart-Vercauteren][2] ciphertext
packing techniques and the [Gentry-Halevi-Smart][3] optimizations. See [this
report][7] for a description of a few of the algorithms using in this library.

Please refer to [CKKS-security.md](CKKS-security.md) for the latest discussion
on the security of the CKKS scheme implementation in HElib.

Since mid-2018 HElib has been under extensive refactoring for *Reliability*,
*Robustness & Serviceability*, *Performance*, and most importantly *Usability*
for researchers and developers working on HE and its uses.

HElib supports an *"assembly language for HE"*, providing low-level routines
(set, add, multiply, shift, etc.), sophisticated automatic noise management,
improved BGV bootstrapping, multi-threading, and also support for Ptxt
(plaintext) objects which mimics the functionality of Ctxt (ciphertext) objects.
The report [Design and implementation of HElib][11] contains additional details.
Also, see [CHANGES.md](CHANGES.md) for more information on the HElib releases.

Full installation instructions and a list of the required dependencies can be
found in [INSTALL.md](INSTALL.md).

For guidance in getting started programming with HElib, take a look at the
example programs and our CKKS tutorials located in the `examples` directory. See
[examples/README.md](examples/README.md).

If you are interested in contributing to HElib, please read our
[Contributing Guidelines](CONTRIBUTING.md).

HElib is written in C++17 and uses the [NTL mathematical library][4].  
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
  [11]: https://eprint.iacr.org/2020/1481    "HElib Design"
  
