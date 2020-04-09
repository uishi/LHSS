#pragma once

#include <cassert>
#include <random>

#include <x86intrin.h>
#include "defines.hpp"

namespace lhss
{

void ConvMpz2Ull(const mpz_class& z, unsigned long long& dest)
{
  mpz_export(&dest, 0, -1, sizeof(dest), 0, 0, z.get_mpz_t());
}

void ConvUll2Mpz(const unsigned long long ull, mpz_class& z)
{
  mpz_import(z.get_mpz_t(), 1, -1, sizeof(ull), 0, 0, &ull);
}

inline std::uint64_t MulUint64Uint64High(const std::uint64_t x, const std::uint64_t y)
{
  return static_cast<std::uint64_t>((static_cast<u128>(x) * y) >> 64);
}

inline void MulUint64Uint64(const std::uint64_t a, const std::uint64_t b, std::uint64_t& low, std::uint64_t& high)
{
  low = _mulx_u64(a, b , reinterpret_cast<unsigned long long int*>(&high));
}

// mul is supposed to be an array of length 2.
inline void MulUint64Uint64(const std::uint64_t a, const std::uint64_t b, std::uint64_t* mul)
{
  mul[0] = _mulx_u64(a, b , reinterpret_cast<unsigned long long int*>(&mul[1]));
}

// mul is supposed to be an array of length 2.
inline void MulUint64Uint64High(const std::uint64_t a, const std::uint64_t b, std::uint64_t& mul)
{
  _mulx_u64(a, b, reinterpret_cast<unsigned long long int*>(&mul));
}

// mul is supposed to be an array of length 2.
inline unsigned char AddUint64Uint64WCarry(const std::uint64_t a, const std::uint64_t b, std::uint64_t& add)
{
  return _addcarry_u64(0, a, b, reinterpret_cast<unsigned long long int*>(&add));
}

inline void AddUint128(const UIntType* a, const UIntType* b, UIntType *res)
{
  unsigned char c = 0;
  c = _addcarry_u64(0, a[0], b[0], reinterpret_cast<unsigned long long int*>(&res[0]));
  _addcarry_u64(c, a[1], b[1], reinterpret_cast<unsigned long long int*>(&res[1]));
}

inline void AddUint128(const UIntType* a, UIntType *res)
{
  unsigned char c = 0;
  c = _addcarry_u64(0, a[0], res[0], reinterpret_cast<unsigned long long int*>(&res[0]));
  _addcarry_u64(c, a[1], res[1], reinterpret_cast<unsigned long long int*>(&res[1]));
}

std::uint32_t BitReverse(std::uint32_t v, std::size_t numbit) 
{
  std::uint32_t r = v; // result
  std::uint32_t s = numbit - 1;
  for (v >>= 1; v; v >>= 1)
  {   
    r <<= 1;
    r |= v & 1; // LSB of current v's 
    s--;
  }
  r <<= s;
  std::uint32_t mask = 0xffffffff;
  std::uint32_t diff = 32 - numbit;
  for (std::uint32_t i = 0; i < diff; ++i)
    mask >>= 1;
  r = r & mask;
  return r;
}

// vector size must be a power of 2.
template <class T>
void PermuteBitRev(const std::vector<T>& src, std::vector<T>& dest)
{
  if (src.size() != dest.size())
    throw std::invalid_argument("vector sizes are different.");
  std::uint32_t len = src.size();
  std::uint32_t numbit = std::log2(len); //equivalently position of highest sinificant bit
  for (std::uint32_t i = 0; i < len; ++i)
    dest[BitReverse(i, numbit)] = src[i];
}

void XGCD(std::int64_t& d, std::int64_t& a, std::int64_t& b, std::int64_t x, std::int64_t y)
{
   std::int64_t u, v, u0, v0, u1, v1, u2, v2, q, r;

   u1 = 1; u2 = 0; u = x; 
   v1 = 0; v2 = 1; v = y;

   while (v != 0)
   {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 = u1 - q * u2;
      v2 = v1 - q * v2;
      u1 = u0;
      v1 = v0;
   }
   d = u;
   a = u1;
   b = v1;
}

// NOTE
// x and mod must be smaller than 2^{63}
void MultInverse(const std::uint64_t x, const std::uint64_t mod, std::uint64_t& inv)
{
  std::uint64_t xx = x;
  if (xx > mod)
    xx = xx % mod;

  std::int64_t a, b, d;
  // d = a*x + b*mod
  XGCD(d, a, b, static_cast<std::int64_t>(xx), static_cast<std::int64_t>(mod));
  // d != 1 means there is no multiplicative inverse modulo mod
  if (d != 1)
    throw std::invalid_argument("Multiplicative inverse is not found.");
  // 1 = a*x mod so, a is equivalent to the inverse
  if (a < 0)
    inv = static_cast<std::uint64_t>(a + mod); 
  else
    inv = static_cast<std::uint64_t>(a);
}

UIntType ExpAndMod(const UIntType base, const UIntType exp, const UIntType modulus)
{
  std::size_t log2_exp_floor = std::log2(exp);
  UIntType nearest_small_pof2 = 1UL << log2_exp_floor;

  // Compute up to base^{log2_exp_floor}
  DoubleUIntType pow = static_cast<DoubleUIntType>(base) * base;
  pow %= modulus;

  std::vector<UIntType> pows{base, static_cast<UIntType>(pow)};
  for (UIntType i = 2; i <= log2_exp_floor; ++i)
  {
    pow *= pow;
    pow %= modulus;
    pows.push_back(pow);
  }
  UIntType exp_acc = nearest_small_pof2;

  UIntType diff, log2_diff_floor;
  // nearest_small_pof2 = 1UL << log2_diff_floor;
  while (exp_acc != exp)
  {
    diff = exp - exp_acc;
    log2_diff_floor = std::log2(diff);

    pow *= pows[log2_diff_floor];
    pow %= modulus;

    nearest_small_pof2 = 1UL << log2_diff_floor;
    exp_acc += nearest_small_pof2;
  }
  return static_cast<UIntType>(pow);
}

// from SEAL
bool IsPrimitiveRoot(const UIntType root, UIntType m, UIntType modulus)
{
  if (root == 0)
    return false;

  // root^{m/2} = -1 mod moulus?.
  UIntType rem = ExpAndMod(root, m/2, modulus);
  return rem == (modulus - 1);
}

// from SEAL
bool FindPrimitiveRoot(const UIntType m, const UIntType prime_modulus, UIntType &root)
{   
  UIntType size_entire_group = prime_modulus - 1;
  UIntType size_quotient_group = size_entire_group / m;
  if (size_entire_group - size_quotient_group * m != 0)
    return false;
  // For randomness
  std::random_device rd; 

  int attempt_counter = 0;
  int attempt_counter_max = 100;
  do
  {
    attempt_counter++;
    // Set destination to be a random number modulo modulus
    root = (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
    root %= prime_modulus;

    // Raise the random number to power the size of the quotient
    // to get rid of irrelevant part
    // std::cout <<  "size_quotient_group=" << size_quotient_group << std::endl;
    root = ExpAndMod(root, size_quotient_group, prime_modulus);
    // std::cerr <<  "root = " << root << std::endl;
  } while (!IsPrimitiveRoot(root, m, prime_modulus) && (attempt_counter < attempt_counter_max));

  return IsPrimitiveRoot(root, m, prime_modulus);
}

void MultModEqualNaive(const std::uint64_t a, const std::uint64_t mod, std::uint64_t& res)
{
  u128 prod = static_cast<u128>(a) * res;
  res = static_cast<std::uint64_t>(prod % mod);  
}

} // namespace