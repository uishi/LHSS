#pragma once

#include <vector>

#include "params.hpp"
#include "z.hpp"

namespace lhss
{
/**
 *  A ring element of \mathbb{Z}[X]/(X^n + 1) for a power-of-two n where its coefficient vector is stored as a member.
 */
struct Poly
{
  std::vector<UIntType> coeffs;

  Poly() : coeffs(Params::n, 0) {}

  // deep copy
  Poly(const std::vector<UIntType>& elements) : coeffs(elements) {}

  UIntType& GetVal(std::size_t i) { return coeffs[i]; }

  const UIntType& GetConstVal(std::size_t i) const { return coeffs[i]; }

  void SetCoeffs(const std::vector<UIntType>& val)
  {
    if (val.size() > coeffs.size())
      throw std::invalid_argument("The input coeffs lenght exceeds the ring dimension.");
    for (std::size_t i = 0; i < val.size(); ++i) 
      coeffs[i] = val[i];
  }

  void SetConstant(const UIntType val) { coeffs[0] = val; }

  void SetDenseSameElem(UIntType val)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      coeffs[i] = val;
  }

  void SetCoeff(const std::size_t index, const UIntType val)
  {
    if (index < Params::n)
      throw std::invalid_argument("index must be < n");

    coeffs[index] = val;
  }

  // TODO delete this once making sure that this is no longer called anywhere.
  void SetMonomial(const std::size_t index, const UIntType val)
  {
    if (index < Params::n)
      throw std::invalid_argument("index must be < n");

    coeffs[index] = val;
  }

  void PackPolyAscend()
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      coeffs[i] = i;
  }

  void PackPolyDescend()
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      coeffs[Params::n - i - 1] = i;
  }

  /**
   * Applys automorphism on a polynomial w.r.t Power Basis.
   * Constant term will remain the same.
   */
  static void AutomorphismNonNTT(const Poly& in, const std::size_t amt, const UIntType mod, Poly& out)
  {
    /*
          a[amt] = a[1];     X -> X^{1 * amt}
      a[2 * amt] = a[2];   X^2 -> X^{2 * amt}
                       ...
      a[k * amt] = a[k];   X^k -> X^{k * amt}
    */
    std::int64_t is_non_zero = 0;
    UIntType res_coeff = 0;
    out.coeffs[0] = in.coeffs[0];
    for (std::size_t i = 1; i < Params::n; ++i)
    {
      std::size_t index_dest = i * amt;
      std::size_t index_dest_mod_n = index_dest & Params::n_minus_one;
      res_coeff = in.coeffs[i];
      // Negate dest coeff if needed
      if ((index_dest >> Params::log2_n) & 1)
      {
        is_non_zero = (in.coeffs[i] != 0);
        res_coeff = (mod - in.coeffs[i]) & static_cast<UIntType>(-is_non_zero);
      }
      out.coeffs[index_dest_mod_n] = res_coeff;
    }
  }

  /**
   * Applys automorphism on a polynomial in NTT domain by an amount amt.
   * SEAL routine.
   */
  static void AutomorphismNTT(const Poly& in, const std::size_t amt, Poly& out)
  {
    UIntType index_dest_revbit_order, index_dest, src_index;
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      index_dest_revbit_order = BitReverse(static_cast<std::uint32_t>(i), Params::log2_n);
      index_dest = amt * (2 * index_dest_revbit_order + 1);
      index_dest &= Params::m_minus_one;
      src_index = BitReverse((index_dest - 1) >> 1, Params::log2_n);
      out.coeffs[i] = in.coeffs[src_index];
    }
  }

  // Use bigger Size?
  static void AddLazy(const Poly& p1, const Poly& p2, Poly& dest)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      dest.coeffs[i] = p1.coeffs[i] + p2.coeffs[i];
  }

  static void AddMod(const Poly& p1, const Poly& p2, UIntType modulus, Poly& p_add)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      ModArith::AddMod(modulus, p1.GetConstVal(i), p2.GetConstVal(i), p_add.coeffs[i]);
  }

  static void AddMod(ModArith& zq, const Poly& p1, const Poly& p2, Poly& p_add)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      zq.AddMod(p1.GetConstVal(i), p2.GetConstVal(i), p_add.coeffs[i]);
  }

  static void AddModEqual(UIntType modulus, const Poly& p, Poly& p_add)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      ModArith::AddModEqual(modulus, p.GetConstVal(i), p_add.coeffs[i]);
  }

  static void AddConstantNTTModEqual(UIntType modulus, UIntType constant, Poly& dest)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      ModArith::AddModEqual(modulus, constant, dest.coeffs[i]);
  }

  static void SubMod(const Poly& p1, const Poly& p2, const ModArith& zq, Poly& p_sub)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      zq.SubMod(p1.GetConstVal(i), p2.GetConstVal(i), p_sub.coeffs[i]);
  }

  // Component-wise mult
  static inline void MulMod(
    const Poly& p1,
    const Poly& p2,
    const ModArith& mod_arith,
    Poly& p_mul
  )
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      mod_arith.MulMod(p1.GetConstVal(i), p2.GetConstVal(i), p_mul.coeffs[i]);
  }

  static inline void MulMod(const ModArith& zq, const Poly& p1, const Poly& p2, Poly& p_mul)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      zq.MulMod(p1.GetConstVal(i), p2.GetConstVal(i), p_mul.coeffs[i]);
  }


  static inline void MulConstantMod(const Poly& p1, const UIntType& c, const ModArith& mod_arith, Poly& p_mult)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      mod_arith.MulMod(p1.GetConstVal(i), c, p_mult.coeffs[i]);
  }

  static inline void MulConstantModEqual(const Poly& p1, const UIntType& c, const ModArith& mod_arith, Poly& p_mult)
  {
    // FIXME Shoup
    for (std::size_t i = 0; i < Params::n; ++i)
      mod_arith.MulMod(p1.GetConstVal(i), c, p_mult.coeffs[i]);
  }

  static void NegateEqual(UIntType modulus, Poly& p)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      ModArith::NegateEqual(modulus, p.coeffs[i]);
  }

  // [0, q) -> [-q/2, q/2)
  static void ReduceToSymmetricEqual(const UIntType& modulus, Poly& p)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      ModArith::ReduceToSymmetricEqual(modulus, p.GetVal(i));
  }

  static std::ostream& ShowSomeCoeff(const Poly& p, const std::size_t beg, const std::size_t end, std::ostream& os)
  {
    os << "[";
    for (std::size_t i = beg; i < end - 1; ++i)
      os << p.coeffs[i] << ",";
 
    os << p.coeffs[end - 1];
    if (end != Params::n)
      os << "...]";      
    else
      os << "]";
  
    return os;
  }

  std::size_t NumNonZeroCoeff()
  {
    std::size_t cnt = 0;
    for (std::size_t i = 0; i < Params::n; ++i)
      if (coeffs[i] != 0) ++cnt;
    
    return cnt;
  }

  static bool IsCoeffsEqual(const Poly& p1, const Poly& p2)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      if (p1.GetConstVal(i) != p2.GetConstVal(i))
        return false;
    
    return true;
  }

  friend std::ostream& operator<<(std::ostream& os, const Poly& p)
  {
    os << "[";
    for (std::size_t i = 0; i < Params::n - 1; ++i)
      os << p.coeffs[i] << ",";
  
    os << p.coeffs[Params::n - 1] << "]";
    return os;
  }
};

struct SignedPoly
{
  std::vector<SIntType> coeffs;

  SignedPoly() : coeffs(Params::n, 0) {}

  SIntType& operator[](std::size_t i) { return coeffs[i]; }

  friend std::ostream& operator<<(std::ostream& os, const SignedPoly& p)
  {
    os << "[";
    for (std::size_t i = 0; i < Params::n - 1; ++i)
      os << p.coeffs[i] << ",";

    os << p.coeffs[Params::n - 1] << "]";
    return os;
  }
};

struct BigPoly
{
  std::vector<mpz_class> coeffs;

  BigPoly() : coeffs(Params::n, 0) {}

  const mpz_class& GetConstVal(std::size_t i) const { return coeffs[i]; }

  static bool IsCoeffsEqual(const BigPoly& p1, const BigPoly& p2)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      if (p1.GetConstVal(i) != p2.GetConstVal(i))
        return false;
    return true;
  }

  static void Add(const BigPoly& p1, const BigPoly& p2, BigPoly& padd)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      padd.coeffs[i] = p1.coeffs[i] + p2.coeffs[i];
  }

  static void AddMod(const BigPoly& p1, const BigPoly& p2, BigPoly& padd)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      padd.coeffs[i] = p1.coeffs[i] + p2.coeffs[i];
      mpz_mod(padd.coeffs[i].get_mpz_t(), padd.coeffs[i].get_mpz_t(), Params::Q.get_mpz_t());
    }
  }

  // For Debugging; Super Slow
  static void Mul(const BigPoly& p1, const BigPoly& p2, BigPoly& pmul)
  {
    mpz_class tmp[2 * Params::n];
    for (std::size_t i = 0; i < 2 * Params::n; ++i)
      tmp[i] = 0;

    mpz_class prod;
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      for (std::size_t j = 0; j < Params::n; ++j)
      {
        prod = p1.coeffs[i] * p2.coeffs[j];
        tmp[i + j] += prod;
        // mpz_mod(tmp[i+j].get_mpz_t(), tmp[i+j].get_mpz_t(), Params::Q.get_mpz_t());
      }
    }

    for (std::size_t i = 0; i < Params::n; ++i)
      pmul.coeffs[i] = tmp[i] - tmp[i + Params::n];
  }

  // For Debugging; Much Slower
  static void MulMod(const BigPoly& p1, const BigPoly& p2, BigPoly& pmul)
  {
    mpz_class tmp[2 * Params::n];
    for (std::size_t i = 0; i < 2 * Params::n; ++i)
      tmp[i] = 0;

    mpz_class prod;
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      for (std::size_t j = 0; j < Params::n; ++j)
      {
        prod = p1.coeffs[i] * p2.coeffs[j];
        mpz_mod(prod.get_mpz_t(), prod.get_mpz_t(), Params::Q.get_mpz_t());
        tmp[i + j] += prod;
        mpz_mod(tmp[i+j].get_mpz_t(), tmp[i+j].get_mpz_t(), Params::Q.get_mpz_t());
      }
    }

    for (std::size_t i = 0; i < Params::n; ++i)
    {
      mpz_class neg_coeff = - tmp[i + Params::n];
      mpz_mod(neg_coeff.get_mpz_t(), neg_coeff.get_mpz_t(), Params::Q.get_mpz_t());
      pmul.coeffs[i] = tmp[i] + neg_coeff;
      mpz_mod(pmul.coeffs[i].get_mpz_t(), pmul.coeffs[i].get_mpz_t(), Params::Q.get_mpz_t());
    }
  }


  // Modulus Reduction to asymmetric interval [0, mod)
  static void ModReduceEqual(const mpz_class& mod, BigPoly& p)
  {
    mpz_class fake;
    for (std::size_t i = 0; i < Params::n; ++i)
      mpz_mod(p.coeffs[i].get_mpz_t(), p.coeffs[i].get_mpz_t(), mod.get_mpz_t());
  }

  static void ScaleAndRound(const BigPoly& p, Poly& p_scaled)
  {
    mpz_class tmp_q, tmp_r, tmp_multby_p;
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      tmp_multby_p = p.GetConstVal(i) * Params::P;

      mpz_tdiv_qr(tmp_q.get_mpz_t(), tmp_r.get_mpz_t(), tmp_multby_p.get_mpz_t(), Params::Q.get_mpz_t());
      
      int sign = mpz_sgn(tmp_q.get_mpz_t());
      int signr = mpz_sgn(tmp_r.get_mpz_t());
      if (Params::z_qi[0].modulus < tmp_q.get_ui())
        throw std::invalid_argument("Sth wrong");

      if (sign > 0)
        p_scaled.coeffs[i] = tmp_q.get_ui();
      else if (sign == 0)
        p_scaled.coeffs[i] = 0;
      else
        p_scaled.coeffs[i] = Params::z_qi[0].modulus - tmp_q.get_ui();
  
      mpz_abs(tmp_r.get_mpz_t(), tmp_r.get_mpz_t());
      // Suppose round(3/7) = round(0.4..), round(4/7) = round(0.5...). Since 7/2 = 3.5, when a numerator >= ceil(7/2) 
      // Negative Case: round(-0.5) = -1.0
      if (tmp_r >= Params::q_half)
      {
        if (signr > 0)
          ++p_scaled.coeffs[i];
        else if (signr < 0)
          --p_scaled.coeffs[i];
      }
      p_scaled.coeffs[i] %= Params::z_qi[0].modulus;
    }
  }

  static void ScaleAndRoundEqual(BigPoly& p)
  {
    mpz_class tmp_r, tmp_multby_p;
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      tmp_multby_p = p.GetConstVal(i) * Params::P;

      mpz_tdiv_qr(p.coeffs[i].get_mpz_t(), tmp_r.get_mpz_t(), tmp_multby_p.get_mpz_t(), Params::Q.get_mpz_t());
      
      int signr = mpz_sgn(tmp_r.get_mpz_t());
      mpz_abs(tmp_r.get_mpz_t(), tmp_r.get_mpz_t());
      // Suppose round(3/7) = round(0.4..), round(4/7) = round(0.5...). Since 7/2 = 3.5, when a numerator >= ceil(7/2) 
      // Negative Case: round(-0.5) = -1.0
      if (mpz_cmp(tmp_r.get_mpz_t(), Params::q_half.get_mpz_t()) >=  0)
      {
        if (signr > 0)
          mpz_add_ui(p.coeffs[i].get_mpz_t(), p.coeffs[i].get_mpz_t(), 1);
        else
          mpz_sub_ui(p.coeffs[i].get_mpz_t(), p.coeffs[i].get_mpz_t(), 1);
      }
      mpz_mod(p.coeffs[i].get_mpz_t(), p.coeffs[i].get_mpz_t(), Params::P.get_mpz_t());
    }
  }

  static void ShowSomeCoeff(const BigPoly& p, const std::size_t beg, const std::size_t end, std::ostream& os)
  {
    os << "[";
    for (std::size_t i = beg; i != end - 1; ++i)
      os << p.coeffs[i] << ",";

    os << p.coeffs[end - 1];
    if (end != Params::n)
      os << "...]";  
    else
      os << "]";
  }

  friend std::ostream& operator<<(std::ostream& os, const BigPoly& p)
  {
    os << "[";
    for (std::size_t i = 0; i < Params::n-1; ++i)
      os << p.GetConstVal(i) << ",";

    os << p.GetConstVal(Params::n - 1) << "]\n";
    return os;
  }
};

struct CRTInt
{
  std::vector<UIntType> rems;

  CRTInt(): rems(Params::num_moduli) {}
};

/* Implementation of R_q \cong R_q1 \times R_q2 \times ... \times R_qk for q = q1 q2 ... qk a product of co-primes */
struct CRTPoly
{
  std::vector<Poly> small_polys;

  CRTPoly() : small_polys(Params::num_moduli) {}

  CRTPoly(const std::size_t k) : small_polys(k) {}

  const Poly& GetConstSmallPoly(std::size_t i) const { return small_polys[i]; }

  void SetConstantPoly(std::vector<UIntType>& crt_constants)
  {
    assert(crt_constants.size() == Params::num_moduli);
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      small_polys[i].coeffs[0] = crt_constants[i];
  }
  
  static void SetConstant(const UIntType val, CRTPoly& p)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
    {
      UIntType v = val;
      if (v >= Params::z_qi[i].modulus)
        v %= Params::z_qi[i].modulus;
  
      p.small_polys[i].coeffs[0] = v;
    }
  }

  static void SetConstantNTT(const UIntType val, CRTPoly& p)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
    {
      UIntType v = val;
      if (v >= Params::z_qi[i].modulus)
        v %= Params::z_qi[i].modulus;

      for (std::size_t j = 0; j < Params::n; ++j)
        p.small_polys[i].coeffs[j] = v;
    }
  }

  static void AutomorphismNTT(const CRTPoly& in, const std::size_t amt, CRTPoly& out)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::AutomorphismNTT(in.GetConstSmallPoly(i), amt, out.small_polys[i]);
  }

  static void AutomorphismNonNTT(const CRTPoly& in, const std::size_t amt, CRTPoly& out)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::AutomorphismNonNTT(in.GetConstSmallPoly(i), amt, Params::z_qi[i].modulus, out.small_polys[i]);
  }


  static void AddMod(const CRTPoly& p1, const CRTPoly& p2, CRTPoly& dest)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::AddMod(p1.GetConstSmallPoly(i), p2.GetConstSmallPoly(i), Params::z_qi[i].modulus, dest.small_polys[i]);
  }

  static void AddMod(const CRTPoly& p, CRTPoly& dest)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::AddModEqual(Params::z_qi[i].modulus, p.GetConstSmallPoly(i), dest.small_polys[i]);
  }

  static void AddModEqual(const CRTPoly& p, CRTPoly& dest)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::AddModEqual(Params::z_qi[i].modulus, p.GetConstSmallPoly(i), dest.small_polys[i]);
  }


  static void AddConstantNTTModEqual(UIntType constant, CRTPoly& dest)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::AddConstantNTTModEqual(Params::z_qi[i].modulus, constant, dest.small_polys[i]);
  }

  static void SubMod(const CRTPoly& p1, const CRTPoly& p2, CRTPoly& dest) 
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::SubMod(p1.GetConstSmallPoly(i), p2.GetConstSmallPoly(i), Params::z_qi[i], dest.small_polys[i]);
  }

  static void MulMod(const CRTPoly& p1, const CRTPoly& p2, CRTPoly& dest)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::MulMod(p1.GetConstSmallPoly(i), p2.GetConstSmallPoly(i), Params::z_qi[i], dest.small_polys[i]);
  }


  static void MulConstantMod(const CRTPoly& p, const UIntType& c, CRTPoly& dest)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::MulConstantMod(p.GetConstSmallPoly(i), c, Params::z_qi[i].modulus, dest.small_polys[i]);
  }

  static void MulConstantMod(const CRTPoly& p, const std::vector<UIntType>& vec_c, CRTPoly& dest) {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::MulConstantMod(p.GetConstSmallPoly(i), vec_c[i], Params::z_qi[i].modulus, dest.small_polys[i]);
  }

  static void MulSmallPolyAndCRTConstant(const Poly& p, const std::vector<UIntType>& c, CRTPoly& p_mult)
  {
    assert(c.size() == Params::num_moduli);
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::MulConstantMod(p, c[i], Params::z_qi[i].modulus, p_mult.small_polys[i]);
  }

  static void NegateEqual(CRTPoly& p)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      Poly::NegateEqual(Params::z_qi[i].modulus, p.small_polys[i]);
  }

  static bool IsPolysEqual(const CRTPoly& p1, const CRTPoly& p2)
  {
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
    {
      if(!Poly::IsCoeffsEqual(p1.GetConstSmallPoly(i), p2.GetConstSmallPoly(i)))
        return false;
    }
    return true;
  }

  friend std::ostream& operator<<(std::ostream& os, const CRTPoly& p)
  {
    os << "[";
    for (std::size_t i = 0; i < Params::num_moduli; ++i)
      os << p.small_polys[i] << "\n";
    os << "]";
    return os;
  }
};

struct PolyConv
{
  static void SymmetricUSmallToISmall(const Poly& upoly, SignedPoly& ipoly)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
      ipoly.coeffs[i] = upoly.coeffs[i];
 }

  // Odd Case  p = 7: {0, 1, 2, 3, -3, -2, -1} \cong {0, 1, 2, 3, 4, 5, 6}
  //         Perform  -7  if val >= 4 = ceil(7/2)
  // Even Case p = 6: {0, 1, 2, -3 -2 -1} \cong {0, 1, 2, 3, 4, 5} 
  //         Perform  -6 if val >= 3 = ceil(6/2)
  static void AsymmetricRepToSymmetricRepEqual(BigPoly& poly)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      if (poly.coeffs[i] >= Params::q_half)
        poly.coeffs[i] -= Params::Q;
    }
  }

  static void AsymmetricRepToSymmetricRepEqual(const mpz_class& mod, BigPoly& poly)
  {
    mpz_class mod_half = mod/2 + 1;
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      if (poly.coeffs[i] >= mod_half)
        poly.coeffs[i] -= mod;
    }
  }

  static void AsymmetricRepToSymmetricRep(const mpz_class& mod, const BigPoly& poly_src, BigPoly& poly_dest)
  {
    mpz_class mod_half = mod/2;
    ++mod_half;
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      poly_dest.coeffs[i] = poly_src.coeffs[i];
      if (poly_dest.coeffs[i] >= mod_half)
        poly_dest.coeffs[i] -= mod;
    }
  }

  static void SymmetricRepToAsymmetricRepEqual(BigPoly& poly)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      if (mpz_sgn(poly.coeffs[i].get_mpz_t()) == -1)
        poly.coeffs[i] += Params::Q;
    }
  }

  static void SymmetricRepToAsymmetricRepEqual(const mpz_class& mod, BigPoly& poly)
  {
    for (std::size_t i = 0; i < Params::n; ++i)
    {
      if (mpz_sgn(poly.coeffs[i].get_mpz_t()) == -1)
        poly.coeffs[i] += mod;
    }
  }


  static void CRTToBigPoly(const CRTPoly& crt_p, BigPoly& p_dest)
  {
    mpz_class rem;
    for (std::size_t j = 0; j < Params::n; ++j)
    {
      for (std::size_t i = 0; i < Params::num_moduli; ++i)
      {
        rem = crt_p.small_polys[i].coeffs[j];
        rem *= Params::crt_q.q_div_qi_invs_mod_qis[i];
        rem *= Params::crt_q.q_div_qis[i];
        p_dest.coeffs[j] += rem;
      }
      // p_dest.coeffs[j] %= Params::Q; // Costly
      mpz_mod(p_dest.coeffs[j].get_mpz_t(), p_dest.coeffs[j].get_mpz_t(), Params::Q.get_mpz_t()); // Costly
    }
  }

  // TODO correct??
  static void BigToCRTPoly(const BigPoly& bp, CRTPoly& p_dest)
  {
    mpz_class fake;
    for (std::size_t j = 0; j < Params::n; ++j)
    {
      for (std::size_t i = 0; i < Params::num_moduli; ++i)
        p_dest.small_polys[i].coeffs[j] = mpz_mod_ui(fake.get_mpz_t(), bp.GetConstVal(j).get_mpz_t(), Params::z_qi[i].modulus);
    }
  }

  // Assuming P = q0
  static void BigToSmallPoly(const BigPoly& bp, Poly& p_dest)
  {
    if (Params::p.size() != 1)
    {
      throw std::invalid_argument("# CRT moduli for p must be 1.");
    }
    mpz_class fake;
    for (std::size_t j = 0; j < Params::n; ++j)
      p_dest.coeffs[j] = mpz_mod_ui(fake.get_mpz_t(), bp.GetConstVal(j).get_mpz_t(), Params::z_qi[0].modulus);
  }
};
} // namespace
