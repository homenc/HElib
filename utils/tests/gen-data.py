#!/usr/bin/env python3

import sys
import argparse
import math
import functools

def sigmoid(x, s):
    return math.exp(x-s)/(1 + math.exp(x-s))

parser = argparse.ArgumentParser()
parser.add_argument("lines", help="total number of lines to generate", type=int)
parser.add_argument("scheme", help="scheme to use", type=str, choices=["BGV", "CKKS"])
parser.add_argument("--columns", help="number of columns to print out", type=int, default=1)
args = parser.parse_args()

if args.scheme == "BGV":
    fn = lambda x: x
elif args.scheme == "CKKS":
    fn = functools.partial(sigmoid, s=args.lines*args.columns/2)

print(args.lines)
for i in range(0, args.lines):
    row = [fn(j) for j in range(i*args.columns, (i+1)*args.columns)]
    print(*row, sep=', ')
