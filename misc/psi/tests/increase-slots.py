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

def parseHeader(line):
    return tuple(int(x) for x in line.split())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input file file", type=str)
    parser.add_argument("nslots", help="what slots length should be", type=int)
    parser.add_argument("--pad-one", help="padding with one (instead of zero)",
        action='store_true', default=False)
    args = parser.parse_args()

    padValue = 1 if args.pad_one is True else 0

    with open(args.infile) as f:
        nrecords, ncols = parseHeader(f.readline())
        print(nrecords, ncols)

        for line in f:
            line = ast.literal_eval(line)
            extendLength = args.nslots - len(line)
            if extendLength > 0:
                line.extend([[padValue]] * extendLength)
            print(repr(line))

if __name__ == "__main__":
    main()
