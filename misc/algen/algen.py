#! /usr/bin/env python3

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

import re
import sys
import argparse
import math
import numth

def parseRange(rangeStr):
  """
  Parses the ranges and numbers given to it.

  Args:
      rangeStr: a string e.g. '2-5,7,10-11'
  Returns:
      A list of values to try out.
  Usage:
      >>> parseRange('2-5,7,10-11')
      [2, 3, 4, 5, 7, 10, 11]
  """
  rangeList = rangeStr.split(',')

  # Group1 regex for range
  # Group2 start number
  # Group3 end number
  # Group4 OR regex for single number
  # Group5 end number

  regex  = re.compile(r'^((\d+)\s*-\s*(\d+))|(\s*(\d+)\s*)$')

  try:
    matches = [ regex.fullmatch(r).groups() for r in rangeList ]
  except AttributeError:
    raise argparse.ArgumentTypeError(\
           "Wrong syntax for range given '%s'. Correct example '2-5,7' " % rangeStr)

  def buildRanges(match):
    matchRange, startRange, endRange, matchSingle, endSingle = match

    if matchRange != None:
      if endRange < startRange:
        raise argparse.ArgumentTypeError(\
          "Range going from high to low '%s'" % matchRange)
      else:
        return range(int(startRange), int(endRange)+1)
    elif matchSingle != None:
      return range(int(endSingle), int(endSingle)+1)
    else:
      raise ValueError(\
              "Something went wrong generating a range with match '%s'"%(match, ) )

  # Transform into a list of range objects
  ranges = list(map(buildRanges, matches))
  # Set Comprehension - guarantee uniqueness
  values = { x for r in ranges for x in r }
  # Convert to sorted list
  values = sorted(list(values))

  return values


def parsePrimeRange(rangeStr):
  """
  Parses the ranges and numbers given to it, but only of primes.

  Args:
      rangeStr: a string e.g. '2-5,7,10-11'
  Returns:
       A list of prime numbers to try out.
  Usage:
      >>> parsePrimeRange('2-5,7,10-11')
      [2, 3, 5, 7, 11]
  """

  p_prime = list(\
              filter(lambda x: numth.factorize(x)[x] == 1, \
              parseRange(rangeStr))\
            )
  if len(p_prime) == 0:
    raise argparse.ArgumentTypeError("No primes found in range given.")

  return p_prime


def algebras(m_maxs):
  """
  Generator for the possible algebras.

  Args:
      m_maxs: Maximum possible m values given p and d as a tuple e.g. (m, p, d).
  Returns:
      A 5-tuple ( m, p, d, phi_m, nSlots ) for the algebra.
  Usage:
      >>> print( '\\n'.join(map(str,algebras( ((24, 5, 2),) ))))
      (2, 5, 1, 1, 1)
      (4, 5, 1, 2, 2)
      (8, 5, 2, 4, 2)
      (3, 5, 2, 2, 1)
      (6, 5, 2, 2, 1)
      (12, 5, 2, 4, 2)
      (24, 5, 2, 8, 4)
  """
  # Generate the divisors for each max_m
  for m_max,p,d in m_maxs:
    factors = numth.factorize(m_max)
    for mFactors in numth.divisorsFactors(factors):
      d_real = d
      m = numth.calcDivisor(mFactors)
      if m == 1:
        continue

      phi_m = numth.phi(mFactors)

      # Correct for order of p in mod m a.k.a. d
      phimFactors = numth.factorize(phi_m)
      for e in sorted(numth.divisors(phimFactors)):
        if p**e % m == 1:
          d_real = e
          break

      nSlots, r = divmod(phi_m, d_real)
      if r != 0:
         raise ArithmeticError("Fractional nslots should not exist.")

      # SolnObj just a 5-tuple ( m, p, d, phi_m, nSlots )
      yield (m,p,d_real,phi_m,nSlots)


def printTable( headers, data, colwidths=None ):
  """
  Print a nice table of possible algebras.

  Args:
      headers: A list/tuple of strings of column headers.
      data: An iterable of list/tuple to print as the row data.
      colwidths: A list/tuple of column widths to be used (forced min. of 8 spaces).
                 Default 12 spaces if None given.
  Returns:
      None.
  Usage:
      >>> printTable(('m','p','d','phi(m)','nSlots'),([(24,5,2,8,4),(12,5,2,4,2)]), [1,1,1,1,1])
         m       p       d     phi(m)  nSlots
            24       5       2       8       4
            12       5       2       4       2
  """
  if colwidths == None:
    colwidths = [12]*len(headers)
  else:
    colwidths = [ width+2 if width > 8 else 8 for width in colwidths ]

  for header, width in zip(headers, colwidths):
    print('{:{align}{width}}'.format(header, align='^', width=width), end='')
  print()

  for row in data:
    for datum, width in zip(row, colwidths):
      print('{:{align}{width}}'.format(datum, align='>', width=width), end='')
    print()


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description="")
  parser.add_argument("-p", type=parsePrimeRange, default="2",
    help="Prime number for plaintext space.")
  parser.add_argument("-d", type=parseRange, default="1",
    help="How many coefficients in a slot (the order of p in ZmStar).")
  parser.add_argument("--self-test", action="store_true",
    help="Run doctest on itself.")
  args = parser.parse_args()

  # Run doctest instead of normal operation
  if args.self_test:
    import doctest
    print(doctest.testmod())
    exit(0)

  # Generator for m_max(s)
  m_maxs = [ (p**d - 1, p, d) for p in args.p for d in args.d ]

  # SolnObj just a 5-tuple ( m, p, d, phi_m, nSlots )
  solns = sorted( algebras(m_maxs) )

  m_max = max(m_maxs)
  minWidths = [ math.ceil(math.log10(t)) for t in m_max ]
  minWidths.append(
    math.ceil(math.log10(numth.phi(numth.factorize(m_max[0])))))
  minWidths.append(minWidths[-1])

  printTable( ('m','p','d','phi(m)','nSlots'), solns, minWidths)
