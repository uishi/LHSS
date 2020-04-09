#pragma once

#include <cmath>
#include <limits>
#include <memory>
// #include <mutex> // for std::call_once and std::once_flag
#include <stdio.h>
#include <vector>

#include "defines.hpp"
#include "crt.hpp"
#include "ntt_params.hpp"
#include "util.hpp"
#include "z.hpp"

namespace lhss
{

namespace const_params
{
constexpr std::size_t kMaxDeg = 65536;
constexpr DoubleUIntType kTwoPow64 = static_cast<DoubleUIntType>(1) << 64;
// Some moduli are borrowed from NFLlib
constexpr UIntType kPrimes[] = {
  4611686018326724609ULL, 4611686018309947393ULL, 4611686018282684417ULL, 4611686018257518593ULL, 4611686018232352769ULL, 4611686018171535361ULL, 4611686018106523649ULL, 4611686018058289153ULL,
  4611686018051997697ULL, 4611686017974403073ULL, 4611686017812922369ULL, 4611686017781465089ULL, 4611686017773076481ULL, 4611686017678704641ULL, 4611686017666121729ULL, 4611686017647247361ULL, 4611686017590624257ULL, 4611686017554972673ULL, 
  4611686012710551553ULL, 4611686012664414209ULL, 4611686012653928449ULL, 4611686012616179713ULL, 4611686012546973697ULL, 4611686012475670529ULL
  };
}

struct Params
{
  // Ring Dimension of Z[X]/(X^n + 1)
  inline static std::size_t n;
  inline static std::size_t log2_n;
  inline static std::size_t n_minus_one;
  inline static std::size_t m;
  inline static std::size_t m_minus_one;

  // Plaintext modulus: a product of co-prime integers 
  // NOTE: p must be larger than maximum magnitude of message for LPR-based HSS (p is not directly related to plaintext modulus)
  inline static std::vector<UIntType> p;
  inline static mpz_class P;

  // Ciphertext modulus q = \prod_{i=1}^{k} q_i
  inline static std::size_t num_moduli;
  inline static std::vector<ModArith> z_qi;
  inline static mpz_class Q;

  // Q/P mod q_i for LPR enc, i.e., the delta.
  inline static std::vector<UIntType> q_div_p_mod_qi;

  // ceil(q/2)
  inline static mpz_class q_half;
  // Q/P for LPR Dec (scaling factor during encryption and homomorphic addition by constant)
  inline static mpz_class q_div_p;

  // Std for Gaussian Distribution
  inline static double sigma;
  // Error bound
  inline static UIntType Berr;

  inline static CRT crt_q;
  inline static CRT crt_p;

  // A primitive m-th (m=2n) root of unity modulo q_i 
  inline static std::vector<UIntType> mth_roots_of_unities;
  // Multiplicative inverse of a primitive m-th root of unity modulo q_i
  inline static std::vector<UIntType> inv_mth_roots_of_unities;

  // For Negative Wrapped-around NTT routines
  inline static std::vector<NTTParams> ntts;

  static std::size_t GetLogQ()
  {
    return mpz_sizeinbase(Q.get_mpz_t(), 2);
  }

  static std::size_t GetLogP()
  {
    return mpz_sizeinbase(P.get_mpz_t(), 2);
  }

  static std::size_t GetN()
  {
    return n;
  }

  /**
   * Assigns proper amount of ciphertext moduli q_i.
   */
  static void AssignModulus()
  {
    z_qi.reserve(num_moduli);
    for (std::size_t i = 0; i < num_moduli; ++i)
    {
      // printf(" -  modulus q%lu = %lu\n", i, const_params::kPrimes[i]);
      z_qi.push_back(ModArith(const_params::kPrimes[i]));
    }
  }

  /**
   * Computes parameters needed for the conversion from Q to P.
   */
  static void ComputeCRTConverter()
  {
    std::vector<UIntType> qi;
    qi.reserve(num_moduli);

    for (auto&& elem: z_qi)
      qi.push_back(elem.modulus);

    crt_q.Set(qi);
    Q = crt_q.q;
  }

  /**
   * Computes Q/P (i.e. bigint) for scaling-up a plaintext message
   */
  static void ComputeQDivP()
  {
    if (p.size() ==  1)
    {
      P = p[0];
    }
    else
    {
      CRT crt_p(p);
      P = crt_p.q;
    }
    // NOTE:
    q_div_p = Q / P;
  }

  /**
   * Computes Q/P \mod q_i for LPR scaling in RNS
   */
  static void ComputeQdivPModQis()
  {
    q_div_p_mod_qi.resize(num_moduli, 1);

    std::vector<UIntType> qis_div_p;

    // for each qi
    for (std::size_t i = 0; i < num_moduli; ++i)
    {
      bool is_divisible = false;
 
      // for each pj
      // NOTE: Since qi is a prime, divisible condition (qi mod pj == 0) is equivalent to (qi == pj)
      for (std::size_t j = 0; j < p.size(); ++j)
      {
        if (z_qi[i].modulus == p[j])
        {
          is_divisible = true;
          break;
        }
      }
      if (!is_divisible)
      {
        qis_div_p.push_back(z_qi[i].modulus);
        q_div_p_mod_qi[i] = 0;
      }
    }

    for (std::size_t i = 0; i < num_moduli; ++i)
    {
      if (q_div_p_mod_qi[i] == 0)
        continue;

      UIntType rem = 1;
      for (auto& q_div_p: qis_div_p)
        z_qi[i].MulModEqual(q_div_p, rem);
    
      q_div_p_mod_qi[i] = rem;
    }

    // ceil(q/2)
    mpz_cdiv_q_ui(q_half.get_mpz_t(), crt_q.q.get_mpz_t(), 2);
  }

  /**
   * Computes powers of 2n-th root of unity mod q_i and those of the inverse of 2n-th root of unity mod q_i for number theoretic transform with negative wrap-around.
   */
  static void ComputeRootsOfUnity()
  {
    // nth_roots_of_unities.reserve(num_moduli);
    // inv_nth_roots_of_unities.reserve(num_moduli);
    mth_roots_of_unities.reserve(num_moduli);
    inv_mth_roots_of_unities.reserve(num_moduli);

    for (std::size_t i = 0; i < num_moduli; ++i)
    {
      UIntType mth_root;
      UIntType inv_mth_root;

      if (!FindPrimitiveRoot(2*n, z_qi[i].modulus, mth_root))
        throw std::invalid_argument("\tA prtimitive 2n-th root of unity was not found.");
    
      MultInverse(mth_root, z_qi[i].modulus, inv_mth_root);
      mth_roots_of_unities.push_back(mth_root);
      inv_mth_roots_of_unities.push_back(inv_mth_root);
    }
  }

  static void ShowParams()
  {
    printf("====== Params List ======\n");
    printf(" -    Ring Dimension = %lu\n", n);
    // printf(" -   Plaintext space = %lu\n", p);
    printf(" - Ciphertext moduli (# moduli = %lu):\n", num_moduli);
    // printf(" - (mod, zeta_n, zeta_n^{-1}, zeta_2n, zeta_2n^{-1}) : ");
    printf(" - (mod, zeta_2n, zeta_2n^{-1}) : ");
    // printf("  {(mod, Barrett Low, Barrett High)} = {");
    for (std::size_t i = 0; i < num_moduli; ++i)
    {
      printf(" \n\t(%lu,", z_qi[i].modulus);
      // printf(" %lu,", nth_roots_of_unities[i]);
      // printf(" %lu,", inv_nth_roots_of_unities[i]);
      printf(" %lu,", mth_roots_of_unities[i]);
      printf(" %lu)", inv_mth_roots_of_unities[i]);
      // printf(" (%lu, %lu, %lu) ",z_qi[i].modulus, z_qi[i].barrett_factor_low, z_qi[i].barrett_factor_high);
    }
    printf("\n");
    gmp_printf(" -         Q: %Zd \n", Q.get_mpz_t());
    printf(" - Q in bits: %lu\n", mpz_sizeinbase(Q.get_mpz_t(), 2));
    printf(" - Sigma = %lf\n", sigma);
    printf("=========================\n");
  }

  // NTT Params: Store powers of (inverse of) a 2n-th root of unity in bit-reveresed order.
  static void ComputeNTTParams()
  {
    ntts.reserve(Params::num_moduli);
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
    {
      UIntType mth_root = Params::mth_roots_of_unities[i];
      UIntType inv_mth_root = Params::inv_mth_roots_of_unities[i];
      ntts.push_back(NTTParams(mth_root, inv_mth_root, Params::z_qi[i], Params::n));
    } // end for i
  }

  /**
   * Assigns Q from user specific RNS moduli
   */
  static void ComputeQ(const std::vector<UIntType>& ctxt_moduli)
  {
    // TODO: co-primality test
    Q = 1;
    z_qi.reserve(ctxt_moduli.size());
    for (auto qi: ctxt_moduli)
    {
      Q *= qi;
      z_qi.push_back(ModArith(qi));
    }
  }

  /**
   * Assigns p from RNS moduli
   */
  static void ComputeP(const std::vector<UIntType>& plain_moduli)
  {
    P = 1;
    p.reserve(plain_moduli.size());
    for (auto plain: plain_moduli)
    {
      p.push_back(plain);
      P *= plain;
    }
  }

  // NOTE: Since qi is a prime, divisible condition (qi mod pj == 0) is equivalent to (qi == pj)
  static bool IsQDivisibleByP(const std::vector<UIntType>& p_moduli)
  {
    std::size_t num_mod_p = p_moduli.size();
    for(std::size_t i = 0; i < num_moduli; ++i)
    {
      for(std::size_t j = 0; j < num_mod_p; ++j)
      {
        if (z_qi[i].modulus == p_moduli[j])
          return true;
      }
    }
    return false;
  }

  static bool IsQDivisibleByP(const UIntType p_modulus)
  {
    for(std::size_t i = 0; i < num_moduli; ++i)
    {
      if (z_qi[i].modulus == p_modulus)
        return true;
    }
    return false;
  }

  static void ComputeBound()
  {
    Berr = std::round(8 * sigma);
  }

  static void Init(
   const std::vector<UIntType>& plain_moduli,
   const std::vector<UIntType>& ctxt_moduli,
   const std::size_t deg
  )
  {
    if (plain_moduli.size() == 0)
      throw std::invalid_argument("No plaintext modulus was found.");

    if (ctxt_moduli.size() == 0)
      throw std::invalid_argument("No ciphertext modulus was found.");
  
    num_moduli = ctxt_moduli.size();
    n = deg;
    log2_n = std::log2(n);
    n_minus_one = n - 1;
    m = n << 1;
    m_minus_one = m - 1;
    
    ComputeQ(ctxt_moduli);
    if (!IsQDivisibleByP(plain_moduli))
      throw std::invalid_argument("Q must be divisible by P.");

    ComputeP(plain_moduli);

    sigma = 3.2;
    ComputeBound();
    // Q
    ComputeCRTConverter();
    ComputeQdivPModQis();
    ComputeQDivP();
    ComputeRootsOfUnity();
    ComputeNTTParams();
  }

  static void Init(const std::size_t nummod, const std::size_t deg, const UIntType plain)
  {
    num_moduli = nummod;
    AssignModulus();

    n = deg;
    log2_n = std::log2(n);
    n_minus_one = n - 1;
    m = n << 1;
    m_minus_one = m - 1;
    
    if (!IsQDivisibleByP(plain))
      throw std::invalid_argument("Q must be divisible by P.");

    p.reserve(1);
    p.push_back(plain);
    P = plain; 
 
    sigma = 3.2;
    ComputeBound();
    // Q
    ComputeCRTConverter();
    ComputeQdivPModQis();
    ComputeQDivP();
    ComputeRootsOfUnity();
    ComputeNTTParams();
  }

  static void Init(const std::size_t nummod, const std::size_t deg, const std::vector<UIntType>& plain)
  {
    num_moduli = nummod;
    AssignModulus();
    if ( ((deg - 1) | deg) != (deg << 1) - 1 || deg <= 0)
      throw std::invalid_argument("Ring dimension must be a power of 2 integer");
  
    n = deg;
    log2_n = std::log2(n);
    n_minus_one = n - 1;
    m = n << 1;
    m_minus_one = m - 1;

    if (!IsQDivisibleByP(plain))
    {
      for (auto& mmod: plain)
        std::cerr << mmod << std::endl;
      throw std::invalid_argument("Q must be divisible by P");
    }
    p = plain;

    sigma = 3.2;
    ComputeBound();
    ComputeCRTConverter();
    ComputeQdivPModQis();
    ComputeQDivP();
    ComputeRootsOfUnity();
    ComputeNTTParams();
  }
};
} // namespace lhss
