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
const std::size_t kMaxDeg = 65536;
constexpr DoubleUIntType kTwoPow64 = static_cast<DoubleUIntType>(1) << 64;
// Some moduli are borrowed from NFLlib
const UIntType kPrimes[] = {
  4611686018326724609ULL, 4611686018309947393ULL, 4611686018282684417ULL, 4611686018257518593ULL, 4611686018232352769ULL, 4611686018171535361ULL, 4611686018106523649ULL, 4611686018058289153ULL,
  4611686018051997697ULL, 4611686017974403073ULL, 4611686017812922369ULL, 4611686017781465089ULL, 4611686017773076481ULL, 4611686017678704641ULL, 4611686017666121729ULL, 4611686017647247361ULL, 4611686017590624257ULL, 4611686017554972673ULL, 
  4611686017529806849ULL, 4611686017517223937ULL, 4611686017496252417ULL, 4611686017489960961ULL, 4611686017439629313ULL, 4611686017429143553ULL, 4611686017401880577ULL, 4611686017376714753ULL, 4611686017290731521ULL, 4611686017246691329ULL,
  4611686017244594177ULL, 4611686017215234049ULL, 4611686017208942593ULL, 4611686017196359681ULL, 4611686017013907457ULL, 4611686016879689729ULL, 4611686016867106817ULL, 4611686016709820417ULL, 4611686016667877377ULL, 4611686016649003009ULL,
  4611686016628031489ULL, 4611686016546242561ULL, 4611686016470745089ULL, 4611686016428802049ULL, 4611686016359596033ULL, 4611686016321847297ULL, 4611686016284098561ULL, 4611686016275709953ULL, 4611686016233766913ULL, 4611686016221184001ULL, 
  4611686016202309633ULL, 4611686016187629569ULL, 4611686016175046657ULL, 4611686016168755201ULL, 4611686016164560897ULL, 4611686016137297921ULL, 4611686016118423553ULL, 4611686016093257729ULL, 4611686015969525761ULL, 4611686015831113729ULL,
  4611686015791267841ULL, 4611686015787073537ULL, 4611686015717867521ULL, 4611686015709478913ULL, 4611686015478792193ULL, 4611686015432654849ULL, 4611686015420071937ULL, 4611686015403294721ULL, 4611686015308922881ULL, 4611686015193579521ULL, 
  4611686015032098817ULL, 4611686015013224449ULL, 4611686015004835841ULL, 4611686014910464001ULL, 4611686014891589633ULL, 4611686014887395329ULL, 4611686014860132353ULL, 4611686014753177601ULL, 4611686014679777281ULL, 4611686014639931393ULL, 
  4611686014583308289ULL, 4611686014545559553ULL, 4611686014532976641ULL, 4611686014476353537ULL, 4611686014459576321ULL, 4611686014407147521ULL, 4611686014402953217ULL, 4611686014388273153ULL, 4611686014381981697ULL, 4611686014312775681ULL, 
  4611686014283415553ULL, 4611686014256152577ULL, 4611686014182752257ULL, 4611686014157586433ULL, 4611686014105157633ULL, 4611686014088380417ULL, 4611686014048534529ULL, 4611686013973037057ULL, 4611686013941579777ULL, 4611686013723475969ULL, 
  4611686013721378817ULL, 4611686013715087361ULL, 4611686013639589889ULL, 4611686013482303489ULL, 4611686013419388929ULL, 4611686013377445889ULL, 4611686013270491137ULL, 4611686013262102529ULL, 4611686013169827841ULL, 4611686013119496193ULL, 
  4611686013106913281ULL, 4611686013092233217ULL, 4611686013085941761ULL, 4611686013067067393ULL, 4611686013010444289ULL, 4611686012949626881ULL, 4611686012947529729ULL, 4611686012909780993ULL, 4611686012840574977ULL, 4611686012834283521ULL, 
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
    // When Q is not divisible by P, one can construct a BFV scheme.
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
   * TODO: Scaled powers of roots of unity for lazy reduction
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
        throw std::invalid_argument("\tA prtimitive 2n-th root of unity was not found!");
    
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

  // NTT Params
  // LN16: Store powers of (inverse of) a 2n-th root of unity in bit-reveresed order.
  // Chen et al. 14: Store powers of (inverse of) a 2n-th (n-th) root of unity.
  static void ComputeNTTParams()
  {
    ntts.reserve(Params::num_moduli);
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
    {
      auto modulus = Params::z_qi[i].modulus;
      // m = 2n
      std::vector<UIntType> power_mth_roots_of_unities(Params::n);
      std::vector<UIntType> power_inv_mth_roots_of_unities(Params::n);
      std::vector<UIntType> power_mth_roots_of_unities_shoup(Params::n);
      std::vector<UIntType> power_inv_mth_roots_of_unities_shoup(Params::n);
    
      UIntType mth_root = Params::mth_roots_of_unities[i];
      UIntType inv_mth_root = Params::inv_mth_roots_of_unities[i];

      UIntType power_mth_root = 1;
      UIntType power_inv_mth_root = 1;
      UIntType power_mth_root_shoup = 1;
      UIntType power_inv_mth_root_shoup = 1;

      power_mth_roots_of_unities[0] = power_mth_root;
      power_inv_mth_roots_of_unities[0] = power_inv_mth_root;
      ModArith::ComputeBarrettConstShoup(modulus, power_mth_root, power_mth_root_shoup);
      ModArith::ComputeBarrettConstShoup(modulus, power_inv_mth_root, power_inv_mth_root_shoup);
      power_mth_roots_of_unities_shoup[0] = power_mth_root_shoup;
      power_inv_mth_roots_of_unities_shoup[0] = power_inv_mth_root_shoup;

      for (std::size_t j = 1; j < Params::n; ++j)
      {
        Params::z_qi[i].MulModEqual(mth_root, power_mth_root);
        Params::z_qi[i].MulModEqual(inv_mth_root, power_inv_mth_root);
        power_mth_roots_of_unities[j] = power_mth_root;
        power_inv_mth_roots_of_unities[j] = power_inv_mth_root;

        ModArith::ComputeBarrettConstShoup(modulus, power_mth_root, power_mth_root_shoup);
        ModArith::ComputeBarrettConstShoup(modulus, power_inv_mth_root, power_inv_mth_root_shoup);
        power_mth_roots_of_unities_shoup[j] = power_mth_root_shoup;
        power_inv_mth_roots_of_unities_shoup[j] = power_inv_mth_root_shoup;

      }

      // std::vector<UIntType> bitrev_power_mth_roots_of_unities(2 * Params::n);
      // std::vector<UIntType> bitrev_power_inv_mth_roots_of_unities(2 * Params::n);
      std::vector<UIntType> bitrev_power_mth_roots_of_unities(Params::n);
      std::vector<UIntType> bitrev_power_inv_mth_roots_of_unities(Params::n);
      std::vector<UIntType> bitrev_power_mth_roots_of_unities_shoup(Params::n);
      std::vector<UIntType> bitrev_power_inv_mth_roots_of_unities_shoup(Params::n);
      
      // Store power of zeta2n in bit-reversed order
      PermuteBitRev(power_mth_roots_of_unities, bitrev_power_mth_roots_of_unities);
      // Store power of inv_zeta2n in bit-reversed order
      PermuteBitRev(power_inv_mth_roots_of_unities, bitrev_power_inv_mth_roots_of_unities);

      PermuteBitRev(power_inv_mth_roots_of_unities_shoup, bitrev_power_inv_mth_roots_of_unities_shoup);
      PermuteBitRev(power_mth_roots_of_unities_shoup, bitrev_power_mth_roots_of_unities_shoup);


      UIntType invn, invn_shoup;
      MultInverse(n, z_qi[i].modulus, invn);
      ModArith::ComputeBarrettConstShoup(modulus, invn, invn_shoup);
      // NTTParams<UIntType> ntt(power_nth_roots_of_unities, power_inv_nth_roots_of_unities, power_mth_roots_of_unities, power_inv_mth_roots_of_unities, bitrev_power_mth_roots_of_unities, bitrev_power_inv_mth_roots_of_unities, invn);
      NTTParams ntt(
       bitrev_power_mth_roots_of_unities,
       bitrev_power_inv_mth_roots_of_unities,
       bitrev_power_mth_roots_of_unities_shoup,
       bitrev_power_inv_mth_roots_of_unities_shoup,
       invn,
       invn_shoup);
      ntts.push_back(std::move(ntt));
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
   * Assign p from RNS moduli
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
