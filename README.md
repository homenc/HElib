HElib
=====

HElib is a software library that implements [homomorphic encryption] [6] (HE). Currently available is an implementation of the [Brakerski-Gentry-Vaikuntanathan] [1] (BGV) scheme, along with many optimizations to make homomorphic evaluation runs faster, focusing mostly on effective use of the [Smart-Vercauteren] [2] ciphertext packing techniques and the [Gentry-Halevi-Smart] [3] optimizations.

At its present state, this library is mostly meant for researchers working on HE and its uses. Also currently it is fairly low-level, and is best thought of as "assembly language for HE". That is, it provides low-level routines (set, add, multiply, shift, etc.), with as much access to optimizations as we can give. Hopefully in time we will be able to provide higher-level routines.

We still do not have an implementation of bootstrapping, hence so far we can only support "leveled" HE, where the parameters must be set large enough to evaluate the desired functionality.

This library is written in C++ and uses the [NTL mathematical library] [4] (version 6.1.0 or higher). It is distributed under the terms of the [GNU General Public License] [5] (GPL).

  [1]: http://eprint.iacr.org/2011/277       "BGV12"
  [2]: http://eprint.iacr.org/2011/133       "SV11"
  [3]: http://eprint.iacr.org/2012/099       "GHS12"
  [4]: http://www.shoup.net/ntl/             "NTL"
  [5]: http://www.gnu.org/licenses/gpl.html  "GPL"
  [6]: http://en.wikipedia.org/wiki/Homomorphic_encryption "Homomorphic encryption"
