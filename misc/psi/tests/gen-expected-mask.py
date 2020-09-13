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
import operator
from enum import Enum
from functools import reduce

def modCoeffsInSlots(listOfLists, p):
    return list([coeff % p for coeff in slot ] for slot in listOfLists)

def modRecord(record, p):
    return list( modCoeffsInSlots(datum, p) for datum in record )

def checkSame(a_iter, b_iter):
    return [ [1] if a == b else [0] for a, b in zip(a_iter, b_iter) ]

def isSameTest(record, querydata):
    return checkSame(record[0], querydata[0])

def andTest(record, querydata):
    mask_a = checkSame(record[0], querydata[0])
    mask_b = checkSame(record[1], querydata[1])
    return [ [1] if a == [1] and b == [1] else [0] for a, b in zip(mask_a, mask_b) ]

def orTest(record, querydata):
    mask_a = checkSame(record[0], querydata[0])
    mask_b = checkSame(record[1], querydata[1])
    return [ [1] if a == [1] or b == [1] else [0] for a, b in zip(mask_a, mask_b) ]

def expandTest(record, querydata):
    mask_a = checkSame(record[0], querydata[0])
    mask_b = checkSame(record[1], querydata[1])
    mask_c = checkSame(record[2], querydata[2])
    return [ [1] if a == [1] or (b == [1] and c == [1]) else [0] for a, b, c in zip(mask_a, mask_b, mask_c) ]

def makeSameSize(x, y):
    maxSz = max(len(x), len(y))
    x.extend([0] * (maxSz - len(x)))
    y.extend([0] * (maxSz - len(y)))
    return (x, y)

def multConstToPtxt(ptxt, const):
    return [ [ const * i for i in slot ] for slot in ptxt ]

def addTwoPtxts(ptxt1, ptxt2):
    return [ list(map(operator.add, *makeSameSize(slot1, slot2))) for slot1, slot2 in zip(ptxt1, ptxt2) ]

def multTwoPtxts(ptxt1, ptxt2):
    return [ list(map(operator.mul, *makeSameSize(slot1, slot2))) for slot1, slot2 in zip(ptxt1, ptxt2) ]

def scoreTest(record, querydata):
    Fs = ( (0, 1, 2, 3),
           (2, 3),
           (0,),
           (1,) )
    mus = (1, 5, 0, 1)
    taus = ( (0, 7, 0, 1),
             (1, 2),
             (1,),
             (1,) )

    masks = [ checkSame(rcol, qcol) for rcol, qcol in zip(record, querydata) ]

    F_scores = []
    for F, tau, mu in zip(Fs, taus, mus):
        mults = (multConstToPtxt(column, t) for column, t in zip((masks[i] for i in F), tau))
        theSum = reduce(addTwoPtxts, mults)
        theSum = addTwoPtxts(theSum, [ [mu] for _ in masks[0] ])
        F_scores.append(theSum)

    return reduce(multTwoPtxts, F_scores)

class TestType:

    funcs = {
               "SAME": isSameTest,
               "AND": andTest,
               "OR": orTest,
               "EXPAND": expandTest,
               "SCORE": scoreTest
            }

    def parse(s):
        for k, v in TestType.funcs.items():
            if s == k:
                return v
        raise ValueError(f"Invalid test option: '{s}'.")

def parseHeader(line):
    return tuple(int(x) for x in line.split())

def readRecord(data, p):
    nrecords, ncols = parseHeader(data.readline())
    for _ in range(nrecords):
        record = [ ast.literal_eval(data.readline()) for _ in range(ncols) ]
        if p > 0:
            record = modRecord(record, p)
        yield record

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("querydata", help="query data file", type=str)
    parser.add_argument("database", help="database file", type=str)
    parser.add_argument("--mod-p", help="mod p", type=int)
    parser.add_argument("--test", help="which test", type=TestType.parse)
    args = parser.parse_args()

    testFn = args.test if args.test else isSameTest
    if not testFn:
        print("Given an invalid test type.", file=sys.stderr)
        sys.exit(1)
    args.mod_p = 0 if args.mod_p is None else args.mod_p

    with open(args.querydata) as f:
        querydata = next(readRecord(f, args.mod_p))

    with open(args.database) as db:
        records = tuple(readRecord(db, args.mod_p))
        print(len(records))
        for record in records:
            print(testFn(record, querydata))

if __name__ == "__main__":
    main()
