# Directory misc/ 

## Introduction

The `misc` directory is essentially a mix of old code that we thought useful to
keep and experimental stuff that may in the future be promoted to a better
suited place in the codebase.

## Contents of `misc/`

The majority of files contained in `misc` are old code and tests, however some
files and directories of note are

### `estimator`
A script that runs the lwe-estimator (https://bitbucket.org/malb/lwe-estimator) and provides approximations of it via linear functions.
Requires `sage` runtime (https://www.sagemath.org/).  

### `params1a` 
An experimental program used to generate various BGV parameters.  This version
only generates `m`s that are products of distinct prime powers.

### `format.sh` 
A script used to format the codebase using `clang-format` greater than or equal to 9. NOTE: This
Ideally the script should be executed from the top directory level `HElib` using `./misc/format.sh
<optional-clang-format>`.

### `aes` 
This directory contains various programs and tests that perform AES
homomorphically using HElib.

### `psi` 
This directory contains experimental code and tests that perform PSI protocols
including database lookup and weighted partial matching algorithms.  View
[README.md](psi/README.md) for more information.

### `algen` 
This directory contains a program `algen.py` that generates BGV parameters.
View [README.md](algen/README.md) for more information.
