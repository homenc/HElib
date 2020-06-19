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

import sys
import argparse
from itertools import islice
import re
import csv

def nSlotsFromParamsOutFile(filename):
    with open(filename) as f:
        pattern = re.compile(r"\s*nslots\s*=\s*(\d+)\s*")
        for line in f:
            match = pattern.fullmatch(line)
            if match:
                return int(match.group(1))
        return None

def printPtxts(ptxts):
    print(len(ptxts))
    for ptxt in ptxts:
        print(ptxt)

def genSlot(csv, scheme):
    fn = int if scheme == "BGV" else float 
    for row in csv:
        yield [ fn(elem) for elem in row ]

def genPtxt(csv, nslots, scheme):
    filler = [0] if scheme == "BGV" else [0.0]
    slots = genSlot(csv, scheme)
    while True:
        ptxt = [ slot for slot in islice(slots, nslots) ]
        if not ptxt:
          break
        ptxt.extend([filler] * (nslots - len(ptxt)))
        yield ptxt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input data file", type=str)
    parser.add_argument("infofile", help="number of slots in a ptxt", type=str)
    parser.add_argument("scheme", help="scheme being used", type=str, 
                        choices=["BGV", "CKKS"])
    args = parser.parse_args()

    # Grab nslots from handy infofile
    nslots = nSlotsFromParamsOutFile(args.infofile)
    if nslots is None:
        sys.exit("Could not find nslots in '%s'." % args.infofile)

    # Read csv file (should just be numbers)
    with open(args.infile) as csv_file:
        next(csv_file) # Skip first line (header)
        csv_reader = csv.reader(csv_file, delimiter=",")
        ptxts = list(genPtxt(csv_reader, nslots, args.scheme))

    printPtxts(ptxts)

if __name__ == '__main__':
    main()
