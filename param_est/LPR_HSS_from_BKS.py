# Code based on the [BKS19, Appendix D.3]
# Finds plain and ciphertext moduli s.t. failure probability per RMS multiplication is at most 2^{-\kappa}.
# User needs to pass the \kappa and B_{max} which is the maximum bound of the message during the computation.

import math

BOUND_SK          = 1
ERROR_STD         = 3.2
BOUND_ERROR       = 8 * ERROR_STD
HAMMING_WEIGHT_SK = 64
INITIAL_RING_DIM  = 8192

def compute_plain_mod(N: int, Bmax: int, hsk: int):
    return  (2 ** (kappa + 2)) * N * Bmax * hsk

def compute_ctxt_mod(p: int, N: int, Bmax: int, Berr: int, hsk: int):
    return (2 ** (kappa + 3)) *  (N ** 2) * Bmax * Berr * (2 * hsk + 1) * p

if __name__ == "__main__" :
    sigma = ERROR_STD
    Berr = BOUND_ERROR
    hsk = HAMMING_WEIGHT_SK
    N = INITIAL_RING_DIM

    Bmax_list  = range(1, 65)
    kappa_list = range(40, 41)

    for e_fail in kappa_list:
      print(" Failure Prob. / RMS Mult. <= 2^-{} ".format(e_fail))
      # exponent of the failure probability of a local rounding
      kappa = e_fail
      for e_bmax in Bmax_list:
        # upper bound of {input, intermediate, output} messages in a circuit
        Bmax = 1 << e_bmax
        print("  Upper bound of the message: 2^{} ".format(e_bmax))

        p = compute_plain_mod(N, Bmax, hsk)

        # minimal q
        q = compute_ctxt_mod(p, N, Bmax, Berr, hsk)
        #qp = compute_ctxt_mod_from_p(p, N, Bmax, Berr, hsk)
        #print(q)

        print("\t (p, logp) = ({}, {})".format(p, math.log(p, 2)))
        print("\t (q, logq) = ({}, {})".format(q, math.log(q, 2)))
