#include <iostream>
#include "lhss.hpp"

using namespace lhss;
using namespace std;

int main ()
{
  // Parameter Setup
  // The maximum magnitude (in bits) of the messages during the secure computation should be provided.
  Setup(32); // message magnitudes are upperbounded by 2^64

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