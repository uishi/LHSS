#pragma once

#include "params.hpp"
#include "poly.hpp"
#include "z.hpp"

namespace lhss
{
#if defined(__clang__)
#define NUMUNROLL 4
#endif
struct SmallNTT
{
  static inline void ForwardButterfly(
    const ModArith& mod_arith,
    const UIntType s,
    const UIntType s_shoup,
    UIntType& x,
    UIntType& y
    )
  {
    UIntType u, v;
    u = x;
    v = y;
    mod_arith.MulModPreCom(s, s_shoup, v);
    mod_arith.AddMod(u, v, x);
    mod_arith.SubMod(u, v, y);
  }

  static inline void BackwardButterfly(
    const ModArith& mod_arith,
    const UIntType phi_inv,
    const UIntType phi_inv_shoup,
    UIntType& x,
    UIntType& y
    )
  {
    UIntType u, v;
    u = x;
    v = y;
    mod_arith.AddMod(u, v, x);
    mod_arith.SubMod(u, v, y);
    mod_arith.MulModPreCom(phi_inv, phi_inv_shoup, y);
  }

  static inline void BackwardButterflyLast(
    const ModArith& mod_arith,
    const UIntType phi_inv,
    const UIntType phi_inv_shoup,
    const UIntType invn,
    const UIntType invn_shoup,
    UIntType& x,
    UIntType& y
    )
  {
    UIntType u, v;
    u = x;
    v = y;
    mod_arith.AddMod(u, v, x);
    mod_arith.SubMod(u, v, y);
    mod_arith.MulModPreCom(invn, invn_shoup, x);
    mod_arith.MulModPreCom(phi_inv, phi_inv_shoup, y);
  }

  /** 
   * Performs forward transformation based on LN16
   * @param[in] param NTT params including powers of a primitive 2n-th root of unity in "bit-reversed ordering"
   + @param[in] zq modular arithmetic over \Z_q 
   + @param[in] p polynomial (w.r.t Power basis) of coefficients in standard ordering
   * @param[out] p polynomial (in NTT domain) of coefficients in bit-reveresed ordering
   */
  static void ApplyNTT(const NTTParams& param, const ModArith& zq, Poly& p)
  {
    std::size_t t = Params::n;
    std::size_t j1, j2;

    for (std::size_t m = 1; m < Params::n; m <<= 1)
    {
      t >>= 1;
      for (std::size_t i = 0; i < m; ++i)
      {
        j1 = 2 * i * t;
        j2 = j1 + t - 1;
        auto s = param.bitrev_pow_zeta2n[m + i];
        auto s_shoup = param.bitrev_pow_zeta2n_shoup[m + i];
  
#if defined(__clang__)
        #pragma unroll NUMUNROLL
#endif
        for (std::size_t j = j1; j <= j2; ++j)
          ForwardButterfly(zq, s, s_shoup, p.coeffs[j], p.coeffs[j+t]);
      }
    }
  }

  /**
   * Performs backward transformation based on LN16
   * @param[in] param NTT params including powers of inverse of a primitive 2n-th root of unity in "bit-reversed ordering"
   + @param[in] zq modular arithmetic over \Z_q 
   * @param[in] p polynomial (in NTT domain) of coefficients in bit-reveresed ordering
   * @param[out] p polynomial (w.r.t Power basis) of coefficients in standard ordering
   */
  static void ApplyInvNTT(const NTTParams& param, const ModArith& zq, Poly& p)
  {
    std::size_t t = 1;
    std::size_t j1, j2, h;
    // auto modulus = zq.modulus;
    // auto two_times_modulus = zq.modulus << 1;
    for (std::size_t m = Params::n; m > 2; m >>= 1)
    {
      j1 = 0;
      h = m >> 1;
      for (std::size_t i = 0; i < h; ++i)
      {
        j2 = j1 + t - 1;
        auto phi_inv = param.bitrev_pow_inv_zeta2n[h + i];
        auto phi_inv_shoup = param.bitrev_pow_inv_zeta2n_shoup[h + i];

#if defined(__clang__)
        #pragma unroll NUMUNROLL
#endif
        for (std::size_t j = j1; j <= j2; ++j)
          BackwardButterfly(zq, phi_inv, phi_inv_shoup, p.coeffs[j], p.coeffs[j + t]);

        j1 += 2*t;
      }
      t <<= 1;
    }

    j1 = 0;
    h  = 1;
#if 1
    for (std::size_t i = 0; i < 1; ++i)
    {
      j2 = j1 + t - 1;
      auto phi_inv_last = param.inv_zeta2n_last_times_invn;
      auto phi_inv_last_shoup = param.inv_zeta2n_last_times_invn_shoup;
      auto invn = param.invn;
      auto invn_shoup = param.invn_shoup;

#if defined(__clang__)
      #pragma unroll NUMUNROLL
#endif
      for (std::size_t j = j1; j <= j2; ++j)
        BackwardButterflyLast(zq, phi_inv_last, phi_inv_last_shoup, invn, invn_shoup, p.coeffs[j], p.coeffs[j + t]);

      j1 += 2*t;
    }
#else
    for (std::size_t i = 0; i < 1; ++i)
    {
      j2 = j1 + t - 1;
      auto phi_inv = param.bitrev_pow_inv_zeta2n[h + i];
      auto phi_inv_shoup = param.bitrev_pow_inv_zeta2n_shoup[h + i];
      for (std::size_t j = j1; j <= j2; ++j)
        BackwardButterfly(zq, phi_inv, phi_inv_shoup, p.coeffs[j], p.coeffs[j + t]);

      j1 += 2*t;
    }
    for (std::size_t j = 0; j < Params::n; ++j)
      zq.MulModPreCom(param.invn, param.invn_shoup, p.coeffs[j]);
#endif
  }
};

struct BigNTT
{
  static void ApplyNTT(CRTPoly& p)
  {
    for (size_t i = 0; i != Params::num_moduli; ++i)
      SmallNTT::ApplyNTT(Params::ntts[i], Params::z_qi[i],  p.small_polys[i]);
  }

  static void ApplyInvNTT(CRTPoly& p)
  {
    for (size_t i = 0; i != Params::num_moduli; ++i)
      SmallNTT::ApplyInvNTT(Params::ntts[i], Params::z_qi[i], p.small_polys[i]);
  }
};
} // namespace lhss