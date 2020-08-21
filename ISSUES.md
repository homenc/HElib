# Known Issues

These issues listed below are known to be current and are expected to be
addressed. 

The function in BinaryArith `fifteenOrLess4Four` is known to not be currently
thread safe. Therefore set NTL's thread pool to 1 to avoid this issue.  The
function `NTL::SetNumThreads` can be used to set the number of threads in the
pool and `NTL::AvailableThreads` gives the current available threads in the
pool.
