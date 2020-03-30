# LHSS (Lattice-based Homomorphic Secret Sharing)

LHSS provides a C++ implementation of Lattice-based homomorphic secret sharing (HSS) based on [BKS19], in which we can build 

* a secure two-party computation with minimal interactions, and

* two server PIR-style applications supporting expressive query and data manipulation.

[BKS19] Elette Boyle and Lisa Kohl and Peter Scholl, "Homomorphic Secret Sharing from Lattices Without FHE", pp. 3 -- 33, vol. 11477, LNCS, Advances in Cryptology -- EUROCRYPT 2019, Springer Berlin Heidelberg, 2019.

# External Dependency

  * gmp

# Running a Demo Program

    mkdir build
    cd build
    cmake ..
    make
    ./demo/sample

# Description

HSS is designed for secrurely evaluating branching program.
The underlying computational model is called a restricted multiplication straightline (RMS) program.

### Data format 

HSS uses two types of data:

 * an encryption of a message multiplied by secret key vector (input value), and

 * a share of a message multiplied by secret key vector (memory value).


## (Im)Possible Operations

* A multiplication on input value and memory value is possible, however, multiplication on memory values is impossible.

* An addition on two shares is almost free.

* An addition on two input values is possible but not recommended since it increases the failure probability of the multiplication.

# Usage

The following shows an example usage of the library.
The library is header-only, which can by used by simply including `lhss.hpp`.

```cpp
#include <iostream>
#include "lhss.hpp"

using namespace lhss;
using namespace std;

int main ()
{
  // Parameter Setup
  // The maximum magnitude (in bits) of the messages during the secure computation should be provided.
  Setup(64); // message magnitudes are upperbounded by 2^64

  // Key-generation internally perfomed in the client instance
  Client c;
  auto public_key = c.GetPk();

  // 2-out-of-2 additive shares of the secret key
  SkShare s0, s1;
  c.ShareSk(s0, s1);

  // Message encryptions 
  uint64_t m1 = 40;
  uint64_t m2 = 50;
  uint64_t m3 = 20;
  HSSCtxt ct1, ct2, ct3;
  HSSEncrypt(public_key, m1, ct1);
  HSSEncrypt(public_key, m2, ct2);
  HSSEncrypt(public_key, m3, ct3);
 
  SkShare t0_out;
  { // the evaluator0 performs the following operations
    SkShare t0;
    Evaluator e0(public_key, s0);
    e0.ConvHSSCtxtToShare(ct1,  t0);
    e0.MultHSSCtxtAndShare(ct2, t0, t0_out);
    e0.ConvHSSCtxtToShare(ct3, t0);
    e0.AddShares(t0, t0_out);
  }

  SkShare t1_out;
  { // the evaluator1 performs the same operations as above, but with a different secret key share
    Evaluator e1(public_key, s1);
    SkShare t1;
    e1.ConvHSSCtxtToShare(ct1,  t1);
    e1.MultHSSCtxtAndShare(ct2, t1, t1_out);
    e1.ConvHSSCtxtToShare(ct3, t1);
    e1.AddShares(t1, t1_out);
  }

  uint64_t res;
  ReconstShares(t0_out, t1_out, res);
  cout << "Result = " << res << endl; // 2020

  return 0;
}
```
# Parameter Selection

The ciphertext and plaintext moduli, and ring dimension are determined based on the followings:

 * Message bound (i.e., the maximum magnitude of a polynomial coefficient)
 
 * The number of multiplications

 * Security level

# Caveat

  PRF implementation is not yet done. This will have a minor impact on the runtime of the multiplications.
