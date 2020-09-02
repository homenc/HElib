# Copyright (C) 2020 The HElib Project Contributors
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

# lwe-estimator.sage
#
# A script that uses the lwe-estimator to estimate the security levels
# for various parameters. This script saves the data in CSV format,
# where the columns are:
#    n, |q|, log2(1/alpha), X=n/log2(1/alpha), Y=security
#
# The script also print to stdout the linear regression line for
# all the X,Y values, and plots the data in PDF files.

import sys
import argparse
import itertools as it

def estimate(hw):
    # Remove all the tests except "usvp" and "dual"
    toskip=("mitm", "arora-gb", "bkw", "dec")

    data = []
    minX = 100000
    maxX = 0

    with open(f"lwe-estimate-{hw}.csv", "w") as csvfile:
        print(f"n,|q|,log2(1/alpha),x=n/log2(1/alpha),security ({hw})", file=csvfile)

        for n in range(2049, 30000, 2048):
            # some data points with sigma*sqrt(2\pi)=8, others with 8*sqrt(n)
            fromLogq = int(round(n/80, 0))
            toLogq = int(round(n/20, 0))
            # Use more q's at the lower range (=higher-security)
            qs = it.chain(range(fromLogq+5,fromLogq+50,10),range(fromLogq,toLogq,50))
            for logQ in qs:
                q = sage.rings.all.RR(2^logQ) +1
                for s2p in [8.0, 8.0*math.sqrt(n)]:
                    alpha = s2p/sage.rings.all.RR(q)

                    # run the estimator to get the complexity of "usvp" and "dual" attacks
                    if isinstance(hw, int):
                        dist=((-1,1), hw) # sparse weight-hw
                    else:
                        dist=(-1,1)       # dense in {-1,0,1}

                    est = estimate_lwe(n, alpha, q, reduction_cost_model=BKZ.sieve, secret_distribution=dist, skip=toskip)

                    # FIXME: can we always assume that "usvp" and "dual" will be there?
                    sec = min(math.log2(est["usvp"]["rop"]), math.log2(est["dual"]["rop"]))
                    if sec < 70 or sec > 270:
                        continue;

                    # record n/log2(1/alpha) and the security
                    loga = -log(alpha)/log(2.0)
                    xx = n / loga
                    minX = min(xx, minX) 
                    maxX = max(xx, maxX) 
                    data.append([xx,sec])

                    # print one line of the CSV file
                    print(f"{n},{logQ},{round(loga,1)},{round(xx,1)},{round(sec,1)}", file=csvfile)

    # Fit the data points to an affine line
    var('a,b')
    model(x) = a*x+b
    line = find_fit(data, model, solution_dict=True)

    # Plot it so you can see that it makes sense
    title = f"{hw}: sec = {round(line[a],1)}X + {round(line[b],1)}"
    p = points(data) + plot(model(a=line[a],b=line[b]), (x,minX,maxX), color='black', title=title)
    p.save(f"lwe-estplot-{hw}.pdf")

    # Also print it to the terminal
    print(title)


def main():
    # Load the estimator (change the location if it is stored elesewhere)
    load("https://bitbucket.org/malb/lwe-estimator/raw/HEAD/estimator.py")
    #load("~/opt/lwe-estimator/estimator.py")

    # Disable information printouts from the estimator
    logging.getLogger("estimator").setLevel(logging.WARNING)

    # Get command-line arguments (if any)
    parser = argparse.ArgumentParser(
        usage='sage %(prog)s [hamming-weight ... ] [--dense]',
        description="""Data is saved in lwe-estimate-<weight>.csv files, and plots  
                       of the data are saved in lwe-estplot-<weight>.pdf files.""")
    parser.add_argument('weights', metavar='hamming-weights', type=int, nargs='*', 
                        help="Hamming weight of keys to estimate. If none specified, defaults to [120,150,...,480] and dense keys.")
    parser.add_argument('--dense', action='store_true', help="estimate (also) dense keys")
    args = parser.parse_args()

    # The default is to run all the tests
    if not args.weights and not args.dense:
        args.weights = list(range(120,500,30))
        args.dense = True # If no HW argument, estimate dense too

    # Get estimates for sparse keys
    if args.weights:
        for hw in args.weights:
            estimate(hw)

    # Get estimates for dense keys
    if args.dense:
        estimate('dense')


if __name__ == "__main__":
    main()
