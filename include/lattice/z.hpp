#pragma once

#include "defines.hpp"
#include "util.hpp"

namespace lhss
{

struct ModArith
{
  UIntType modulus;
  UIntType barrett_const[2]{0, 0}; // \floor(2^{128}/p)

  ModArith(UIntType mod) : modulus(mod) { Init(); }

  // FIXME remove gmp
  void Init()
  {
    mpz_class two_pow_128;
    mpz_ui_pow_ui(two_pow_128.get_mpz_t(), 2, 128);
    mpz_class mpz_mod(modulus);
    mpz_class quotient;
    mpz_fdiv_q(quotient.get_mpz_t(), two_pow_128.get_mpz_t(), mpz_mod.get_mpz_t());
    mpz_class quotient_high;
    mpz_fdiv_q_2exp(quotient_high.get_mpz_t(), quotient.get_mpz_t(), 64);
    barrett_const[0] = mpz_get_ui(quotient.get_mpz_t());
    barrett_const[1] = mpz_get_ui(quotient_high.get_mpz_t());
  }

  static void ComputeBarrettConstShoup(
   const UIntType modulus,
   const UIntType c,
   UIntType& c_barr
   )
  {
    u128 mult = (static_cast<u128>(1) << 64) * c;
    c_barr = mult / modulus;
  }

  // propagating the carries
  inline void ModReduce128(const UIntType* a, UIntType& res) const
  {
    UIntType tmp1, tmp2[2], tmp3, carry;
    MulUint64Uint64High(a[0], barrett_const[0], carry);
    MulUint64Uint64(a[0], barrett_const[1], tmp2);
    tmp3 = tmp2[1] + AddUint64Uint64WCarry(tmp2[0], carry, tmp1);
    MulUint64Uint64(a[1], barrett_const[0], tmp2);
    carry = tmp2[1] + AddUint64Uint64WCarry(tmp1, tmp2[0], tmp1);
    tmp1 = a[1] * barrett_const[1] + tmp3 + carry;
    res = a[0] - tmp1 * modulus;
    res -= modulus & (static_cast<UIntType>(-static_cast<SIntType>(res >= modulus)));
  }

  void AddLazy(const UIntType v1, const UIntType v2, UIntType& dest) const
  {
    dest = v1 + v2;
  }

  void AddMod(const UIntType v1, const UIntType v2, UIntType& dest) const
  {
    dest = v1 + v2;
    if (dest > modulus)
      dest -= modulus;
  }

  void AddModEqual(const UIntType v, UIntType& dest)
  {
    dest = dest + v;
    if (dest > modulus)
      dest -= modulus;
  }

  static void AddMod(const UIntType mod, const UIntType v1, const UIntType v2, UIntType& v_add)
  {
    v_add = v1 + v2;
    if (v_add > mod)
      v_add -= mod;
  }

  static void AddModEqual(const UIntType mod, const UIntType v, UIntType& dest)
  {
    dest = dest + v;
    if (dest > mod)
      dest -= mod;
  }

  void SubMod(UIntType v1, UIntType v2, UIntType& dest) const
  {
    dest = v1 - v2;
    if (v1 < v2)
      dest += modulus;
  }

  static void SubMod(UIntType mod, UIntType v1, UIntType v2, UIntType& dest) 
  {
    dest = v1 - v2;
    if (v1 < v2)
      dest += mod;
  }

  inline void MulLazy(const UIntType v1, const UIntType v2, UIntType* dest) const
  {
    MulUint64Uint64(v1, v2, dest[0], dest[1]);
  }

  inline void MulMod(const UIntType v1, const UIntType v2, UIntType& dest) const 
  {
    UIntType mul[2];
    MulLazy(v1, v2, mul);
    ModReduce128(mul, dest);
  }

  // c \in [0, 4p)
  // v \in [0, 2p)
  // Intermediate value \in [0, 4p)
  // dest \in [0, 2p)
  inline void MulModPreComLazy(
   const UIntType c,
   const UIntType c_precon,
   const UIntType v,
   UIntType& dest
  ) const 
  {
    UIntType approx_quot = MulUint64Uint64High(v, c_precon); 
    approx_quot *= modulus;
    UIntType mult = c * v;
    dest = mult - approx_quot;
  }

  // c \in [0, 4p)
  // v \in [0, 2p) \rightarrow  [0, 2p)
  // Intermediate value \in [0, 4p)
  inline void MulModPreComLazy(
   const UIntType c,
   const UIntType c_precon,
   UIntType& v
  ) const 
  {
    UIntType approx_quot = MulUint64Uint64High(v, c_precon); 
    approx_quot *= modulus;
    UIntType mult = c * v;
    v = mult - approx_quot;
  }

  inline void MulModPreCom(
   const UIntType c,
   const UIntType c_precon,
   UIntType& v
  ) const 
  {
    MulModPreComLazy(c, c_precon, v);
    v -= modulus & (static_cast<UIntType>(-static_cast<SIntType>(v >= modulus)));
  }


  inline void MulModPreCom(
   const UIntType c,
   const UIntType c_precon,
   const UIntType v,
   UIntType& res
  ) const 
  {
    res = v;
    MulModPreComLazy(c, c_precon, res);
    res -= modulus & (static_cast<UIntType>(-static_cast<SIntType>(res >= modulus)));
  }

  inline void MulModEqual(const UIntType v, UIntType& dest) const
  {
    UIntType mul[2];
    MulUint64Uint64(v, dest, mul[0], mul[1]);
    ModReduce128(mul, dest);
  }

  // static void MulModLazy(UIntType v1, UIntType v2, DoubleUIntType& dest)
  // {
  //   dest = static_cast<DoubleUIntType>(v1) * v2;
  // }

  static void NegateEqual(UIntType mod, UIntType& dest)
  {
    dest = mod - dest;
  }

  static void ReduceToSymmetricEqual(const UIntType& mod, UIntType& dest)
  {
    if (dest > ((mod - 1) >> 1))
      dest -= mod;
  }
};

} // namespace lhss
