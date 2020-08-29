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
import ast
from functools import reduce
from itertools import islice
import re

def nSlotsFromParamsOutFile(filename):
    with open(filename) as f:
        pattern = re.compile(r"\s*nslots\s*=\s*(\d+)\s*")
        for line in f:
            match = pattern.fullmatch(line)
            if match:
                return int(match.group(1))
        return None

def parsePtxtGetSlots(iterable):
    for ptxt in iterable:
        for slot in ast.literal_eval(ptxt):
            yield slot

def printoutNums(ptxts, nelements):
    for slot in islice(parsePtxtGetSlots(ptxts), nelements):
        print(*slot, sep=", ")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input data file", type=str)
    parser.add_argument("infofile", help="file containing algebra info",
                        type=str)
    parser.add_argument("scheme", help="scheme being used", type=str,
                        choices=["BGV", "CKKS"])
    parser.add_argument("--nelements", help="total number of elements",
                        type=int)
    args = parser.parse_args()

    nslots = nSlotsFromParamsOutFile(args.infofile)
    if nslots is None:
        sys.exit("Could not find nslots in '%s'." % args.infofile)

    with open(args.infile) as data:
        header = data.readline()
        # Set if nelements is set otherwise calculate from info file.
        nelements = args.nelements if args.nelements is not None \
            else reduce((lambda x, y: x * y), map(int, header.split())) * nslots
        print(nelements)
        printoutNums(data, nelements)

if __name__ == '__main__':
    main()
