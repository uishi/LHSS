#pragma once

#include <cstring>

#include "ctxt.hpp"
#include "secretkey.hpp"
#include "rlwe_ops_sym.hpp"
#include "rlwe_ops.hpp"

namespace lhss
{
class Evaluator
{
 public:
  using SkShare = SecretKey;
  using SkPtr = std::shared_ptr<SkShare>;
  using PkPtr = std::shared_ptr<Ciphertext>;

  Evaluator(
   PkPtr& pk,
   SkShare& sk_share, 
   std::string&& prf_key = "",
   std::size_t party_index = 0,
   bool is_load_opt = false
  ) :
   pk_(*pk),
   sk_share_(sk_share),
   prf_key_(prf_key),
   party_index_(party_index),
   is_load_opt_(is_load_opt)
   {
    if (is_load_opt_)
    { 
      if (!IsSkOptShareVaild())
        throw std::invalid_argument("Sk Share is invalid.");
    }
   }

  Evaluator(Ciphertext& pk, SkPtr sk_share, 
   std::string&& prf_key = "", std::size_t party_index = 0,
   bool is_load_opt = false) 
   : pk_(pk), sk_share_(*sk_share),
     prf_key_(prf_key), party_index_(party_index),
     is_load_opt_(is_load_opt)
   {
    if (is_load_opt_)
      if (!IsSkOptShareVaild())
        throw std::invalid_argument("Sk Share is invalid.");
   }


  /**
   * Checks if coefficients of secret key is less than 2 in magnitude.
   * Eace coefficient is expected to be in {-1, 0, 1}.
   **/
  bool IsSkOptShareVaild()
  {
    for(std::size_t i = 0; i < Params::num_moduli; ++i)
    {
      for(std::size_t j = 0; j < Params::n; ++j)
      {
        if (sk_share_.GetConstB().GetConstSmallPoly(i).GetConstVal(j) >= 2)
          return false;
      }
    }
    return true;
  }

  /**
   * Loads an input value into memory by decrypting HSS ciphertext by sk share
   * @param[in]  ct a HSS ciphertext that encrypts x
   * @param[out]  t a share of (x, x*s)
   */
  void ConvHSSCtxtToShare(const HSSCtxt& c, SkShare& t) const
  {
    BigPoly bp[2];
    HSSDotProd(sk_share_, c, t.elements);
    PolyConv::CRTToBigPoly(t.elements[0], bp[0]);
    PolyConv::CRTToBigPoly(t.elements[1], bp[1]);
    ScaleAndRound(bp);
    Lift(bp, t);
  }

  void MultHSSCtxtAndShare(const HSSCtxt& c, const SkShare& s, SkShare& t_out) const
  {
    BigPoly bp[2];

    HSSDotProd(s, c, t_out.elements);
    PolyConv::CRTToBigPoly(t_out.elements[0], bp[0]);
    PolyConv::CRTToBigPoly(t_out.elements[1], bp[1]);

    ScaleAndRound(bp);
    Lift(bp, t_out);
  }

  /**
   * Applys an automorphism on a share
   * @param[in] in an input share (m(X), m(X)s(X))
   * @param[in] amt amount for automorphism
   * @param[out] out an output share of the form (m(X^(amt)), m(X^(amt)) * s(X^(amt))
   */
  void ApplyAurotmophismNTT(const SkShare& in, const std::size_t amt, SkShare& out) const
  {
    CRTPoly::AutomorphismNTT(in.GetConstB(), amt, out.GetB());
    CRTPoly::AutomorphismNTT(in.GetConstA(), amt, out.GetA());
  }

#if 0
  /**
   * Performs the key-switching on a share. 
   * This is a bit more efficient than the key-switching on the ciphertext.
   * @param[in] in an input share (m(X^i), m(X^i)s(X^i))
   * @param[in] evk hint for the key-swithcing
   * @param[out] out an output share of the form (m(X^i), m(X^i) * s(X))
   */  
  void KeySwitchOnShare() const
  {
  }
#endif

  void AddShares(const SkShare& s1, const SkShare& s2, SkShare& s_out) const
  {
    SkShare::AddMod(s1, s2, s_out);
  }

  void AddShares(const SkShare& s1, SkShare& s_out) const
  {
    SkShare::AddMod(s1, s_out);
  }

 private:
  // TODO implement via AES-NI
  void PRF();

  Ciphertext pk_;
  SkShare sk_share_;
  std::string prf_key_;
  std::size_t party_index_ = 0;
  bool is_load_opt_ = false;
};
} // namespace lhss