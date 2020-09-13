# Program for generating algebras for BGV.

## Contents in directory

* algen.py -- the program to generate BGV algebras.
* numth.py -- utility functions.

Developed with Python 3.6.9, but may work on lower.

## Optional dependency

algen works faster if it uses Linux `factor` utility for factorizing integers.

On mac, you can add it with home brew. It gets installed as `gfactor` algen
detects it.

## How to use

Run
```
./algen -h 
```
for options.

But essentially, you are deciding which prime numbers `p` and how many
polynomial coefficients you want in a slot `d`.

The nice thing is that you can provide ranges (inclusive that is [a,b]) as 
```
./algen -p 2,5 -d 1,4
```

## Future plan

Rewrite this in C++. 
