# Key Generation and Management Pipeline

## Introduction
This subdirectory provides additional utilities for use within a key management system. This code should be considered experimental.

## Difference to other utils
This folder contains a plaintext generator `gen-data.cpp`, as well as `key-gen.cpp` an alternative key generation to `create-context.cpp`. `key-gen.cpp` differs from `create-context.cpp` in the following respects:
-  `key-gen` only generates keys for BGV without bootstrapping, resulting in a smaller executable.
-  `key-gen` writes only the secret key polynomial to the secret key file, rather than printing both the associated public key and the secret key.
-  `key-gen` creates two public keys - `Enc.pk`, a public key for encryption, and `Eval.pk`, a public key which can be used for homomorphic function evaluation. `Enc.pk` is significantly smaller than `Eval.pk`.

## Installation 
Installation is identical to the larger utils directory, and proceeds via running  CMake and then make from the build directory. This will create five executables in `build/bin`, namely `create-context`, `decrypt`, `encrypt`, `gen-data`, and `key-gen`.

## Running the utilities
All utilities have a help method by passing the `-h` flag. In the `key-gen` folder, three parameter sets have been included for demonstration, `small_`, `medium_` and `large_params.txt`. We describe an example pipeline using the `key-gen` utilities from within the key-gen folder, with the small parameters as an example. This will create all files inside the key-gen folder.
1. Generate some test data using the `gen-data` program
```
./../build/bin/gen-data small_params.txt -o small_data -n 5
```
The optional `-o` flag defines the output file name, while the optional `-n` flag determines the number of plaintexts to generate. Running this command will create a file `small_data.json` with five plaintexts in the key-gen directory.
2. Generate keys using `key-gen`
```
./../build/bin/key-gen small_params.txt -s -o small
```
The optional `-s` flag means that only the secret key polynomial will be written, resulting in a smaller file. The optional `-o` flag specifies the prefix of the generated keys: in this case, `small`. This will create three key files: `small.sk`, a secret key, `smallEnc.pk`, a public key for encryption, and `smallEval.pk`, a public key for homomorphic evaluation.

All optional flags that can be passed to `create-context` can also be passed to `key-gen` -- to see a list, run the help method `-h`.
3. Encrypt data
```
./../build/bin/encrypt smallEnc.pk small_data.json -o small.ctxt
```
The optional `-o` argument determines the file name of the resulting ciphertext file. Running this command will encrypt `small_data.json` with public key `smallEnc.pk` to give 5 ciphertexts in the file `small.ctxt`.
4. Decrypt the data
```
./../build/bin/decrypt small.sk small.ctxt -o small_result.json -s
```
This decrypts `small.ctxt` with the secret key `small.sk` and puts the resulting plaintext into `small_result.json`. The `-s` flag indicates that only the secret key polynomial is stored in `small.sk`, and must be used if keys were generated using the `-s` flag with `key-gen`.
