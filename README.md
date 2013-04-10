HElib
=====

HElib is a software library that implements homomorphic encryption (HE). Currently available is an implementation of the [Brakerski-Gentry-Vaikuntanathan] [1] (BGV) scheme, along with many optimizations to make homomorphic evaluation runs faster, focusing mostly on effective use of the [Smart-Vercauteren] [2] ciphertext packing techniques and the [Gentry-Halevi-Smart] [3] optimizations.

This library is written in C++ and uses the [NTL mathematical library] [4]. It is distributed under the terms of the [GNU General Public License] [5] (GPL).

  [1]: http://eprint.iacr.org/2011/277       "BGV12"
  [2]: http://eprint.iacr.org/2011/133       "SV11"
  [3]: http://eprint.iacr.org/2012/099       "GHS12"
  [4]: http://www.shoup.net/ntl/             "NTL"
  [5]: http://www.gnu.org/licenses/gpl.html  "GPL"
