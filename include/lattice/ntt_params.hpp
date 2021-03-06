#pragma once

#include "defines.hpp"
#include "util.hpp"
#include "z.hpp"

namespace lhss
{

struct NTTParams
{
  using UIntVec = std::vector<UIntType>;

  UIntVec bitrev_pow_zeta2n;
  UIntVec bitrev_pow_inv_zeta2n;
  UIntVec bitrev_pow_zeta2n_shoup;
  UIntVec bitrev_pow_inv_zeta2n_shoup;

  UIntType invn = 0;
  UIntType invn_shoup = 0;
#if 1
  // for performing 1/n at the last butterfly together with multiplication by W^-1
  UIntType inv_zeta2n_last_times_invn = 0;
  UIntType inv_zeta2n_last_times_invn_shoup = 0;
#endif

  NTTParams(){}

  NTTParams(
   const UIntType zeta2n,
   const UIntType inv_zeta2n,
   const ModArith& mod_arith,
   const size_t n
  )
  {
    Init(zeta2n, inv_zeta2n, mod_arith, n);
  }

  NTTParams(
   UIntVec& bitrev_ms,
   UIntVec& bitrev_ims,
   UIntVec& bitrev_ms_shoup,
   UIntVec& bitrev_ims_shoup,
   UIntType invn_mod,
   UIntType invn_shoup_mod,
   UIntType inv_zeta2n_last_times_invn_mod,
   UIntType inv_zeta2n_last_times_invn_shoup_mod
  )
   : bitrev_pow_zeta2n(std::move(bitrev_ms)),
     bitrev_pow_inv_zeta2n(std::move(bitrev_ims)),
     bitrev_pow_zeta2n_shoup(std::move(bitrev_ms_shoup)),
     bitrev_pow_inv_zeta2n_shoup(std::move(bitrev_ims_shoup)),
     invn(invn_mod),
     invn_shoup(invn_shoup_mod),
     inv_zeta2n_last_times_invn(inv_zeta2n_last_times_invn_mod),
     inv_zeta2n_last_times_invn_shoup(inv_zeta2n_last_times_invn_shoup_mod)
  {}

  NTTParams(const NTTParams&& move)
  { 
    bitrev_pow_zeta2n = std::move(move.bitrev_pow_zeta2n);
    bitrev_pow_inv_zeta2n = std::move(move.bitrev_pow_inv_zeta2n);
    bitrev_pow_zeta2n_shoup = std::move(move.bitrev_pow_zeta2n_shoup);
    bitrev_pow_inv_zeta2n_shoup = std::move(move.bitrev_pow_inv_zeta2n_shoup);
    invn = move.invn;
    invn_shoup = move.invn_shoup;
    inv_zeta2n_last_times_invn = move.inv_zeta2n_last_times_invn;
    inv_zeta2n_last_times_invn_shoup = move.inv_zeta2n_last_times_invn_shoup;
  }

  // NTTParams(const NTTParams& copy)
  // { 
  //   bitrev_pow_zeta2n = copy.bitrev_pow_zeta2n;
  //   bitrev_pow_inv_zeta2n = copy.bitrev_pow_inv_zeta2n;
  //   bitrev_pow_zeta2n_shoup = copy.bitrev_pow_zeta2n_shoup;
  //   bitrev_pow_inv_zeta2n_shoup = copy.bitrev_pow_inv_zeta2n_shoup;
  //   invn = copy.invn;
  //   invn_shoup = copy.invn_shoup;
  // }

  void Init(
   const UIntType zeta2n,
   const UIntType inv_zeta2n,
   const ModArith& mod_arith,
   const size_t n
  ) 
  {
    UIntType modulus = mod_arith.modulus;

    bitrev_pow_zeta2n.resize(n);
    bitrev_pow_inv_zeta2n.resize(n);
    bitrev_pow_zeta2n_shoup.resize(n);
    bitrev_pow_inv_zeta2n_shoup.resize(n);

    MultInverse(n, modulus, invn);

    UIntType pow_zeta2n = 1;
    UIntType pow_inv_zeta2n = 1;
    UIntType pow_zeta2n_shoup = 1;
    UIntType pow_inv_zeta2n_shoup = 1;

    ModArith::ComputeBarrettConstShoup(modulus, pow_zeta2n, pow_zeta2n_shoup);
    ModArith::ComputeBarrettConstShoup(modulus, pow_inv_zeta2n, pow_inv_zeta2n_shoup);
    UIntVec pow_zeta2ns(n);
    UIntVec pow_inv_zeta2ns(n);
    UIntVec pow_zeta2ns_shoup(n);
    UIntVec pow_inv_zeta2ns_shoup(n);

    pow_zeta2ns[0]           = pow_zeta2n;
    pow_inv_zeta2ns[0]       = pow_inv_zeta2n;
    pow_zeta2ns_shoup[0]     = pow_zeta2n_shoup;
    pow_inv_zeta2ns_shoup[0] = pow_inv_zeta2n_shoup;

    for (std::size_t j = 1; j < n; ++j)
    {
      mod_arith.MulModEqual(zeta2n, pow_zeta2n);
      mod_arith.MulModEqual(inv_zeta2n, pow_inv_zeta2n);
      ModArith::ComputeBarrettConstShoup(modulus, pow_zeta2n, pow_zeta2n_shoup);
      ModArith::ComputeBarrettConstShoup(modulus, pow_inv_zeta2n, pow_inv_zeta2n_shoup);

      pow_zeta2ns[j]           = pow_zeta2n;
      pow_inv_zeta2ns[j]       = pow_inv_zeta2n;
      pow_zeta2ns_shoup[j]     = pow_zeta2n_shoup;
      pow_inv_zeta2ns_shoup[j] = pow_inv_zeta2n_shoup;
    }

    PermuteBitRev(pow_zeta2ns, bitrev_pow_zeta2n);
    PermuteBitRev(pow_inv_zeta2ns, bitrev_pow_inv_zeta2n);
    PermuteBitRev(pow_zeta2ns_shoup, bitrev_pow_zeta2n_shoup);
    PermuteBitRev(pow_inv_zeta2ns_shoup, bitrev_pow_inv_zeta2n_shoup);

    MultInverse(n, mod_arith.modulus, invn);
    ModArith::ComputeBarrettConstShoup(modulus, invn, invn_shoup);
      
    mod_arith.MulMod(bitrev_pow_inv_zeta2n[1], invn, inv_zeta2n_last_times_invn);
    ModArith::ComputeBarrettConstShoup(modulus, inv_zeta2n_last_times_invn, inv_zeta2n_last_times_invn_shoup);
  }
};
} // namespace lhss