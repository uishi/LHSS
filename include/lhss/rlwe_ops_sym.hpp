#pragma once

#include <cassert>

#include "ctxt.hpp"
#include "sampler.hpp"
#include "secretkey.hpp"

namespace lhss
{
/**
 * Performs distributed decrypton on a ciphertext
 * @param[in] sk a shared of secret key (in R^2_Q) in NTT form
 * @param[in] c a ciphertext (in R^2_Q) in NTT form
 * @param[out] dest_plain secret share of output plaintext in Non-NTT form
 */
void HSSDDecryption(const SecretKey& sk, const Ciphertext& c, CRTPoly& dest_plain)
{
  CRTPoly tmp;
  CRTPoly::MulMod(c.GetConstB(), sk.GetConstB(), tmp);
  CRTPoly::MulMod(c.GetConstA(), sk.GetConstA(), dest_plain);
  CRTPoly::AddModEqual(tmp, dest_plain);
  BigNTT::ApplyInvNTT(dest_plain);
}

/**
 * Performs optimal distributed decrypton
 * @param[in] sk a shared of secret key (in R^2_Q) in NTT form
 * @param[in] c a ciphertext (in R^2_Q) in NTT form
 * @param[out] dest_plain secret share of output plaintext in Non-NTT form
 */
void HSSDDecryptionOpt(const SecretKey& sk, const Ciphertext& c, CRTPoly& dest_plain)
{
  CRTPoly::MulMod(c.GetConstA(), sk.GetConstA(), dest_plain);
  if (sk.GetConstB().GetConstSmallPoly(0).GetConstVal(0) == 1)
    CRTPoly::AddModEqual(c.GetConstB(), dest_plain);
  // if not 1, sk.b should be 0 and therefore we do not have to do anything more.
  BigNTT::ApplyInvNTT(dest_plain);
}

void HSSDDecryption(const SecretKey& sk, const HSSCiphertext& c, CRTPoly t[2])
{
  CRTPoly t0, t1;
  CRTPoly::MulMod(c.GetConstB0(), sk.GetConstB(), t0);
  CRTPoly::MulMod(c.GetConstA0(), sk.GetConstA(), t[0]);
  CRTPoly::AddModEqual(t0, t[0]);
  BigNTT::ApplyInvNTT(t[0]);

  CRTPoly::MulMod(c.GetConstB1(), sk.GetConstB(), t1);
  CRTPoly::MulMod(c.GetConstA1(), sk.GetConstA(), t[1]);
  CRTPoly::AddModEqual(t1, t[1]);
  BigNTT::ApplyInvNTT(t[1]);
}

void DotProd(const SecretKey& sk, const Ciphertext& c, CRTPoly& res)
{
  CRTPoly t;
  CRTPoly::MulMod(c.GetConstB(), sk.GetConstB(), t);
  CRTPoly::MulMod(c.GetConstA(), sk.GetConstA(), res);
  CRTPoly::AddModEqual(t, res);
  BigNTT::ApplyInvNTT(res);
}

void HSSDotProd(const SecretKey& sk, const HSSCtxt& c, CRTPoly t[2])
{
  DotProd(sk, c.enc_m, t[0]);
  DotProd(sk, c.enc_m_times_s, t[1]);
}

void HSSDDecryptionOpt(const SecretKey& sk, const HSSCiphertext& c, CRTPoly t[2])
{
  // Multiplication on A
  CRTPoly::MulMod(c.GetConstA0(), sk.GetConstA(), t[0]);
  CRTPoly::MulMod(c.GetConstA1(), sk.GetConstA(), t[1]);
  if (sk.GetConstB().GetConstSmallPoly(0).GetConstVal(0) == 1) 
  {
    CRTPoly::AddModEqual(c.GetConstB0(), t[0]);
    CRTPoly::AddModEqual(c.GetConstB1(), t[1]);
  } // if not 1, sk.b should be 0 and therefore we do not have to do anything more.
  BigNTT::ApplyInvNTT(t[1]);
  BigNTT::ApplyInvNTT(t[0]);
}

/**
 * Performs LPR encryption of zero under the secret key 
 * 
 * @param[in] m plaintext to be encrypted
 * @param[out] ct LPR ciphetext
 */
void EncryptZero(const SecretKey& sk, Ciphertext& ct)
{
  // a
  CRTSampler::SampleUniformModQ(ct.GetA());
  // a*s
  CRTPoly as;
  CRTPoly::MulMod(ct.GetConstA(), sk.GetConstA(), as);
  // e
  CRTPoly e;
  CRTSampler::SamplePolyFromNormalDist(e);
  BigNTT::ApplyNTT(e);
  // e - a*s
  CRTPoly::SubMod(e, as, ct.GetB());
}

/**
 * Performs LPR encryption under the secret key
 * @param[in] m a plaintext to be encrypted
 * @param[out] ct a ciphertext that encrypts m
 */
void Encrypt(const SecretKey& sk, const Poly& m, Ciphertext& ct)
{
  CRTPoly m_times_delta;
  CRTPoly::MulSmallPolyAndCRTConstant(m, Params::q_div_p_mod_qi, m_times_delta);
  BigNTT::ApplyNTT(m_times_delta);

  EncryptZero(sk, ct);
  CRTPoly::AddModEqual(m_times_delta, ct.GetB());
}

/**
 * Performs LPR encryption under the secret key
 * @param[in] m a plaintext to be encrypted
 * @param[out] ct a ciphertext that encrypts m
 */
void EncryptMTimesSk(const SecretKey& sk, const Poly& m, Ciphertext& ct)
{
  CRTPoly m_times_delta;
  CRTPoly::MulSmallPolyAndCRTConstant(m, Params::q_div_p_mod_qi, m_times_delta);
  BigNTT::ApplyNTT(m_times_delta);

  EncryptZero(sk, ct);
  CRTPoly::AddModEqual(m_times_delta, ct.GetA());
}

/**
 * Performs LPR encryption under the secret key
 * @param[in] m a plaintext to be encrypted
 * @param[out] ct a ciphertext that encrypts m
 */
void Encrypt(const SecretKey& sk, const std::uint64_t& m, Ciphertext& ct)
{
  Poly plain;
  plain.SetCoeffs({m});
  CRTPoly m_times_delta;
  CRTPoly::MulSmallPolyAndCRTConstant(plain, Params::q_div_p_mod_qi, m_times_delta);
  BigNTT::ApplyNTT(m_times_delta);

  EncryptZero(sk, ct);
  CRTPoly::AddModEqual(m_times_delta, ct.GetB());
}

/**
 * Performs HSS LPR encryption under the secret key
 * @param[in] m a plaintext to be encrypted
 * @param[out] ct a ciphertext that encrypts m
 */
void HSSEncrypt(const SecretKey& sk, const std::uint64_t& m, HSSCtxt& ct)
{
  Poly plain;
  plain.SetCoeffs({m});
  Encrypt(sk, plain, ct.enc_m);
  EncryptMTimesSk(sk, plain, ct.enc_m_times_s);
}

//FIXME
void CircularEncrypt(const SecretKey& sk, const SecretKey& plain_sk, HSSCiphertext& ct)
{
  // deep copy
  SecretKey org_sk = plain_sk;
  CRTPoly org_s = org_sk.GetA();
  // Convert it to non-NTT rep.
  // BigNTT::ApplyInvNTT(org_s);
  // Multiply Delta
  CRTPoly scaled_plain;
  CRTPoly::MulConstantMod(org_s, Params::q_div_p_mod_qi, scaled_plain);
  // BigNTT::ApplyNTT(scaled_plain);

  CRTPoly e0, e1;
  CRTSampler::SamplePolyFromNormalDist(e0);
  CRTSampler::SamplePolyFromNormalDist(e1);
  BigNTT::ApplyNTT(e0);
  BigNTT::ApplyNTT(e1);  
  // a0, a1
  CRTSampler::SampleUniformModQ(ct.GetA0());
  CRTSampler::SampleUniformModQ(ct.GetA1());

  // a0*s, a1*s
  CRTPoly a0s, a1s;
  CRTPoly::MulMod(ct.GetConstA0(), sk.GetConstA(), a0s);
  CRTPoly::MulMod(ct.GetConstA1(), sk.GetConstA(), a1s);

  CRTPoly::SubMod(e0, a0s, ct.GetB0());
  CRTPoly::SubMod(e1, a1s, ct.GetB1()); 
  // e - a*s + (q/p)m
  CRTPoly::AddModEqual(scaled_plain, ct.GetB0());
  CRTPoly::AddModEqual(scaled_plain, ct.GetA1());
}

// input should be NTTed vals
void ReconstShare(const SecretKey& s0, const SecretKey& s1, BigPoly recons[2])
{
  SecretKey t0 = s0;
  SecretKey t1 = s1;
  SecretKey debug_m;
  SecretKey::AddMod(t0, t1, debug_m);

  BigNTT::ApplyInvNTT(debug_m.GetB());
  BigNTT::ApplyInvNTT(debug_m.GetA());

  PolyConv::CRTToBigPoly(debug_m.GetConstB(), recons[0]);
  PolyConv::CRTToBigPoly(debug_m.GetConstA(), recons[1]);
  // "Output  (x)"
  BigPoly::ModReduceEqual(Params::P, recons[0]);
  PolyConv::AsymmetricRepToSymmetricRepEqual(recons[0]);
  // BigPoly::ShowSomeCoeff(debug_message[0], 0, 10, cout);

  // "Output(x*s)"
  BigPoly::ModReduceEqual(Params::P, recons[1]);
  PolyConv::AsymmetricRepToSymmetricRepEqual(Params::P, recons[1]);
  // BigPoly::ShowSomeCoeff(debug_message[1], 0, 10, cout);
}

// Only the constant term of the message is obtained.
void ReconstShares(const SecretKey& s0, const SecretKey& s1, std::uint64_t& res)
{
  BigPoly recons;
  SecretKey debug_m;
  SecretKey::AddMod(s0, s1, debug_m);

  BigNTT::ApplyInvNTT(debug_m.GetB());
  PolyConv::CRTToBigPoly(debug_m.GetConstB(), recons);
  res = mpz_mod_ui(recons.coeffs[0].get_mpz_t(), recons.coeffs[0].get_mpz_t(), Params::z_qi[0].modulus);
}
} // namespace lhss