# Contributing to the Homomorphic Encryption Library HElib

Adding new features, improving documentation, fixing bugs, writing new 
tests, designing and coding new examples or writing tutorials are all examples 
of helpful contributions. 

A contribution to HElib can be initiated through GitHub pull request (PR).
HElib is written in C++14 and uses `clang-format` for the formatting of the code.
Requiremens and installation instructions can be found in [INSTALL.md](INSTALL.md).
When making code contributions to HElib, we ask that you follow the `C++14` 
coding standard and format your code using the [clang format](.clang-format) 
style file included in this distribution. Please provide unit/regression 
tests that are relevant to your code contribution. 

This project uses Developer Certificate of Origin [DCO](https://developercertificate.org/). 
Be sure to sign off your commits using the `-s` flag or adding `Signed-off-By: Name<Email>` in the commit message.

### Example commit message
```bash
git commit -s -m 'Informative commit message'
```

### Unit/Regression tests

Whether you are contributing a new feature, updating or bug fixing the code elsewhere, you need to ensure that your HElib build passes the regression test suite. 

HElib test suite uses the [Google Test Framework](https://github.com/google/googletest). Additional information on HElib's test suite can be found in [TESTS.md](TESTS.md). Please remember to provide unit/regression tests that are relevant to your code contribution.

Once all the tests have passed, and you are satisfied with your contribution, open a pull request into the `master` branch from **your fork of the repository** to request adding your contributions into the main code base.

