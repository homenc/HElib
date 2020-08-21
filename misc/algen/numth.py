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

import math
import itertools
import shutil
import subprocess
from collections import Counter
from collections import defaultdict
from functools import reduce

# Global var of factor path
factorPath = shutil.which('factor') \
             or shutil.which('gfactor')

# Optimisations or different sieves could be used
def nextPrime(upto):
  """A simple generator implementation of the Sieve of Eratosthenes.

  Args:
      upto: An integer for the upper bound of the sieve.
  Returns:
      A generator which yields an integer the next prime.
  Usage:
      >>> list(nextPrime(10))
      [2, 3, 5, 7]
  """
  numbers = [True]*(upto+1)
  numbers[0] = False
  numbers[1] = False

  p = 2
  yield p # Return 2 straight away

  while p**2 <= upto:
    for i in range(p, math.floor(upto/p)+1):
      numbers[p*i] = False
    for i, n in enumerate(numbers[p+1:], p+1):
      if n == True:
        p = i
        break
    yield p

  # Other primes not encoutered during sieving.
  for i, n in enumerate(numbers[p+1:], p+1):
    if n == True:
      yield i

def twoandOdds(upto):
  # Using range seems fatser that itertools.count
  return itertools.chain((2,), range(3, upto, 2))


def factorize(n):
  """Factorize an integer.

  Args:
      n: An integer to factorize.
  Returns:
      A dict of prime factors (keys) and their exponents (values).
  Usage:
      >>> factorize(1546)
      defaultdict(<class 'int'>, {2: 1, 773: 1})
  """

  if n < 1:
    raise ValueError("No prime factors for %s" % n)

  factors = defaultdict(int)
  if n == 1:
    return factors # Empty dict - nothing to iterate through

  # Factor path global - run once
  if factorPath != None:
    res = subprocess.run(["gfactor", str(n)], stdout=subprocess.PIPE, universal_newlines=True)
    # ignore first numer - the dividend
    factors = map( int, res.stdout.split()[1:] )
    return defaultdict(int, Counter(factors))
  else:
    # This uses far less memory than using nextPrime
    for p in twoandOdds(math.ceil(n/2)+1):
      while n % p == 0 and n != 1:
        n /= p
        factors[p] += 1
      if n == 1:
        break

    if not any(factors):
      factors[n] = 1 # Must be prime.

  return factors


def phi(factors):
  """
  The Euler totient calculated by given prime factors.

  Args:
      factors: Dictionary of prime factors (keys) and their exponents (values).
  Returns:
      An integer. The value of the Euler totient.
  Usage:
      >>> phi({2: 2, 5: 2})
      40
      >>> phi({2: 2, 5: 2, 7: 1})
      240
  """

  phiN = 1

  for factor, count in factors.items():
    if count == 0:
      continue
    phiN = phiN*(factor-1)*factor**(count-1)

  return phiN


def calcDivisor(factors):
  """
  Return the number from its prime factors.

  Args:
      factors: Dictionary of prime factors (keys) and their exponents (values).
  Returns:
      An integer. The number made of the prime factors.
  Usage:
      >>> calcDivisor({5: 3, 7: 2})
      6125
      >>> calcDivisor({2: 3, 5: 2, 7: 1})
      1400
  """

  acc = 1
  for f,e in factors.items():
    acc *= f**e
  return acc


def divisorsFactors(factors):
  """
  A generator from the divisors generated from the prime factors given.

  Args:
      factors: Dictionary of prime factors (keys) and their exponents (values).
  Returns:
      A divisor as a dict of prime factors (keys) and their exponents (values).
  Usage:
      >>> print('\\n'.join(map(str, divisorsFactors({2:2, 5:2}))))
      {2: 0, 5: 0}
      {2: 1, 5: 0}
      {2: 2, 5: 0}
      {2: 0, 5: 1}
      {2: 1, 5: 1}
      {2: 2, 5: 1}
      {2: 0, 5: 2}
      {2: 1, 5: 2}
      {2: 2, 5: 2}
  """

  counter = [0]*len(factors)

  acc = reduce(lambda x,y: x*(y+1), factors.values(), 1)

  for _ in range(acc):
    q = 1 # To kick-off the counting
    yield dict( zip( factors.keys(), counter )) # up front to give zeroth value
    for i,e in enumerate(factors.values()):
      if q == 1:
        q, counter[i] = divmod(counter[i]+1, e+1)
      else:
        break


def divisors(factors):
  """
  A generator that returns the divisors from the prime factors of the dividend.

  Args:
      factors: Dictionary of prime factors (keys) and their exponents (values).
  Returns:
      An integer divisor.
  Usage:
      >>> list(divisors({2: 2, 5: 2}))
      [1, 2, 4, 5, 10, 20, 25, 50, 100]
  """
  for d in divisorsFactors(factors):
    yield calcDivisor(d)


if __name__ == "__main__":
  import doctest
  print(doctest.testmod())
