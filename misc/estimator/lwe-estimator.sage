# bootstrapping-lwe-estimator.sage
#
# A script that uses the lwe-estimator to estimate the security levels
# for various bootstrapping parameters. This file prints the results
# to stdout in CSV format, where the columns are:
#    m, n, |q|, X=n/log2(1/alpha), Y=security
#
# The script also print *to stderr* the linear regression line for
# all the X,Y values.
#
# The bootsrapping parameters are specified by the ms, ns, and qs
# lists below, where n[i]=phi(m[i]) and q[i] is the bit-length of
# the modulus in the key-switching matrices, namely
#   |q| = log2(ctxtPrimes)+log2(specialPrimes).

import sys
import itertools as it

print("Usage: sage", sys.argv[0], "[HammingWeight]\n")
print("    if HammingWeight is not supplied, will use [120,150,180,...,480]")
print("    and dense keys. If HammingWeight='dense' only dense keys are used.\n")
print("    Data is saved in lwe-estimate-<weight>.csv files, and plots")
print("    of the data are saved in lwe-estplot-<weight>.pdf files.")

if len(sys.argv) > 1:
  if sys.argv[1].isnumeric() and int(sys.argv[1])>0:
    weights = [int(sys.argv[1])]
  else:
    weights = ['dense']
else:
  weights=[120]
#  weights=range(120,500,30)

# We have the two estimates sec = 3.8*x -15 and sec=2.53*x +19,
# for dense and sparse(120), respectively, where
# x = n/log_2(1/alpha) = n/log_2(sigma*sqrt(2 pi)/q).
#     n/x = log_2(1/alpha) = log_2(q) -log_2(sigma*sqrt(2 pi))
#     q = 2^{n/x} * sigma*sqrt(2 pi)
#
# To set the range of |q| values on which to run the estimator, we aim
# at security roughly between 80 and 256, so we roughly parameters that
# give x values between 20 and 80.

# Load the estimator (change the location if it is stored elesewhere)
# load("https://bitbucket.org/malb/lwe-estimator/raw/HEAD/estimator.py")
load("~/opt/lwe-estimator/estimator.py")

# Disable information printouts from the estimator
logging.getLogger("estimator").setLevel(logging.WARNING)

# Remove all the tests except "usvp" and "dual"
toskip=("mitm", "arora-gb", "bkw", "dec")

def estimate_hw(hw):
  data = []
  minX = 100000
  maxX = 0
  csvfile = open("lwe-estimate-"+str(hw)+".csv", "w")
  print("n,|q|,log2(1/alpha),x=n/log2(1/alpha),security ("+str(hw)+")", file=csvfile)
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
          dist=((-1,1),hw) # sparse weight-hw
        else:
          dist=(-1,1)      # dense in {-1,0,1}

        est = estimate_lwe(n, alpha, q, reduction_cost_model=BKZ.sieve,secret_distribution=dist,skip=toskip)

        # FIXME: can we always assume that "usvp" and "dual" will be there?
        sec = min(math.log2(est["usvp"]["rop"]), math.log2(est["dual"]["rop"]))
        if sec < 70 or sec > 270:
          continue;

        # record n/log2(1/alpha) and the security
        loga = -log(alpha)/log(2.0)
        xx = n / loga
        if xx < minX:
          minX = xx
        if xx > maxX:
          maxX = xx
        data.append([xx,sec])

        # print one line of the CSV file
        print(str(n)+','+str(logQ)+','+str(round(loga,1))+','+str(round(xx,1))+','+str(round(sec,1)), file=csvfile)

  csvfile.close()

  # Fit the data points to an affine line
  var('a,b')
  model(x)=a*x+b
  line=find_fit(data,model,solution_dict=True)

  # Plot it so you can see that it makes sense
  title = str(hw)+": sec = "+str(round(line[a],1))+"X + "+str(round(line[b],1))
  p=points(data)+plot(model(a=line[a],b=line[b]),(x,minX,maxX),color='black',title=title)
  p.save('lwe-estplot-'+str(hw)+'.pdf')

  # Also print it to the terminal
  print(title)


for hw in weights:
  estimate_hw(hw)

if len(sys.argv) < 2: # no HW argument, estimate dense as well
  estimate_hw('dense')
