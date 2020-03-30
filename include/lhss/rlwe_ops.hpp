#pragma once

#include <cassert>

#include "ctxt.hpp"
#include "sampler.hpp"
#include "secretkey.hpp"

namespace lhss
{
using PkPtr = std::shared_ptr<Ciphertext>;

void MergeDDecPlain(const Ciphertext& c, const std::vector<CRTPoly>& ddec_plains, BigPoly& m)
{
  CRTPoly sum = c.GetConstB();
  BigNTT::ApplyInvNTT(sum);  
  
  for (auto& pl: ddec_plains)
    CRTPoly::AddModEqual(pl, sum);
  
  PolyConv::CRTToBigPoly(sum, m);
  PolyConv::AsymmetricRepToSymmetricRepEqual(m);
  BigPoly::ScaleAndRoundEqual(m);
}

/**
 * Performs LPR encryption of zero under the public key
 * @param[in] m a plaintext to be encrypted
 * @param[out] ct a ciphertext that encrypts m
 */
void EncryptZero(const PkPtr& pk, Ciphertext& ct)
{
  // e0
  CRTSampler::SamplePolyFromNormalDist(ct.GetB());
  // e1
  CRTSampler::SamplePolyFromNormalDist(ct.GetA());
    
  CRTPoly v;
#if SAMPLE_OPT
  CRTSampler::SampleUniformTernary(v);
#else
  CRTSampler::SamplePolyFromNormalDist(v);
#endif

  BigNTT::ApplyNTT(ct.GetB());
  BigNTT::ApplyNTT(ct.GetA());
  BigNTT::ApplyNTT(v);

  // pk.b * v and pk.a * v
  CRTPoly tmp_b, tmp_a;
  CRTPoly::MulMod(pk->GetConstB(), v, tmp_b);
  CRTPoly::MulMod(pk->GetConstA(), v, tmp_a);

  CRTPoly::AddModEqual(tmp_b, ct.GetB());
  CRTPoly::AddModEqual(tmp_a, ct.GetA());
}

/**
 * Performs LPR encryption under the public key
 * @param[in] m a plaintext to be encrypted
 * @param[out] ct a ciphertext that encrypts m
 */
void Encrypt(const PkPtr& pk, const Poly& m, Ciphertext& ct)
{
  CRTPoly m_times_delta;
  CRTPoly::MulSmallPolyAndCRTConstant(m, Params::q_div_p_mod_qi, m_times_delta);
  BigNTT::ApplyNTT(m_times_delta);

  EncryptZero(pk, ct);
  CRTPoly::AddModEqual(m_times_delta, ct.GetB());
}


/**
 * Performs LPR encryption under the public key
 * @param[in] m a plaintext to be encrypted
 * @param[out] ct a ciphertext that encrypts m
 */
void EncryptMTimesSk(const PkPtr& pk, const Poly& m, Ciphertext& ct)
{
  CRTPoly m_times_delta;
  CRTPoly::MulSmallPolyAndCRTConstant(m, Params::q_div_p_mod_qi, m_times_delta);
  BigNTT::ApplyNTT(m_times_delta);

  EncryptZero(pk, ct);
  CRTPoly::AddModEqual(m_times_delta, ct.GetA());
}

/**
 * Performs LPR-based HSS encryption under the public key by generating (Enc(m), Enc(m*s)) \in (R_q)^4.
 * @param[in] pk a public key
 * @param[in] m a plaintext to be encrypted (the infinity norm of the message should be small enough for efficiency)
 * @param[out] lct a HSS ciphertext that encrypts m (and m*s)
 */
// void HSSEncrypt(const PkPtr& pk, const Poly& m, HSSCtxt& hss_ct)
void HSSEncrypt(const PkPtr& pk, const Poly& m, HSSCtxt& hss_ct)
{
  Encrypt(pk, m, hss_ct.enc_m);
  EncryptMTimesSk(pk, m, hss_ct.enc_m_times_s);
}

/**
 * Performs LPR-based HSS encryption under the public key by generating (Enc(m), Enc(m*s)) \in (R_q)^4.
 * @param[in] pk a public key
 * @param[in] m a plaintext to be encrypted (the infinity norm of the message should be small enough for efficiency)
 * @param[out] lct a HSS ciphertext that encrypts m (and m*s)
 */
// void HSSEncrypt(const PkPtr& pk, const std::uint64_t m, HSSCtxt& hss_ct)
void HSSEncrypt(const PkPtr& pk, const std::uint64_t m, HSSCtxt& hss_ct)
{
  Poly plain;
  plain.SetCoeffs({m});
  Encrypt(pk, plain, hss_ct.enc_m);
  EncryptMTimesSk(pk, plain, hss_ct.enc_m_times_s);
}

void ScaleAndRound(BigPoly* bp)
{
  PolyConv::SymmetricRepToAsymmetricRepEqual(bp[0]);
  PolyConv::SymmetricRepToAsymmetricRepEqual(bp[1]);
  BigPoly::ScaleAndRoundEqual(bp[0]);
  BigPoly::ScaleAndRoundEqual(bp[1]);
  BigPoly::ModReduceEqual(Params::P, bp[0]);
  BigPoly::ModReduceEqual(Params::P, bp[1]);
}

void Lift(BigPoly* bp, SecretKey& t)
{
  PolyConv::AsymmetricRepToSymmetricRepEqual(Params::P, bp[0]);
  PolyConv::AsymmetricRepToSymmetricRepEqual(Params::P, bp[1]);
  PolyConv::BigToCRTPoly(bp[0], t.GetB());
  PolyConv::BigToCRTPoly(bp[1], t.GetA());
  BigNTT::ApplyNTT(t.GetB());
  BigNTT::ApplyNTT(t.GetA());
}
} // namespace lhss