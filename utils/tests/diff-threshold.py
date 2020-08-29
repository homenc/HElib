#!/usr/bin/env python3

# Copyright (C) 2020 IBM Corp.
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

import argparse
import sys
import ast
import math

def diff_float(na, nb, threshold):
    for a, b in zip(na, nb):
        if not math.isclose(a, b, abs_tol=threshold):
            raise ValueError(f"Difference {a - b} between {a} and {b} "
                             f"exceeds threshold {threshold}.")

def makeSameSize(a, b, max_length):
    lenA, lenB = len(a), len(b)
    if lenA > max_length or lenB > max_length:
        raise ValueError(f"Size of slots for {a}({lenA}) {b}({lenB}) "
                         f"> {max_length}.")

    if lenA == lenB:
        return a, b
    else:
        maxSz = max(lenA, lenB)
        a += [0] * (maxSz - lenA)
        b += [0] * (maxSz - lenB)
        return (a, b)

def parseCorrectly(la, lb, decrypt):
    error_msg = "Type mismatch. {0}({1}) and {2}({3}) type do not match."
    if decrypt:
        for a, b in zip(la, lb):
            a, b = ast.literal_eval(a), ast.literal_eval(b)
            if type(a) is not type(b):
                raise TypeError(error_msg.format(a, type(a), b, type(b)))
            yield a, b
    else:
        for a, b in zip(la, lb):
            a = [[ float(i) for i in a.split(",") ]]
            b = [[ float(i) for i in b.split(",") ]]
            if type(a) is not type(b):
                raise TypeError(error_msg.format(a, type(a), b, type(b)))
            yield a, b

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("firstfile", help="first data file", type=str)
    parser.add_argument("secondfile", help="second data file", type=str)
    parser.add_argument("--decrypt", help="diff decrypt format (instead of decode)",
                        action='store_true')
    parser.add_argument("--threshold", help="error threshold [default=0.001]",
                        type=float, default=0.001)
    args = parser.parse_args()

    with open(args.firstfile, 'r') as f1, open(args.secondfile, 'r') as f2:
        l1, l2 = list(f1), list(f2)

    if len(l1) != len(l2):
        sys.exit(f"Different number of lines. "
                 f"First contains {len(l1)} second contains {len(l2)}.")

    if l1[0] != l2[0]:
        sys.exit(f"File headers differ. {l1[0]} {l2[0]}.")

    try:
        for a, b in parseCorrectly(l1[1:], l2[1:], args.decrypt):
            for sa, sb in zip(a, b):
                sa, sb = makeSameSize(sa, sb, 2)
                diff_float(sa, sb, args.threshold)
    except (TypeError, ValueError) as e:
        sys.exit(str(e))

if __name__ == "__main__":
    main()
