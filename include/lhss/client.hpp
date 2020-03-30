#pragma once

#include "ctxt.hpp"
#include "sampler.hpp"
#include "secretkey.hpp"

namespace lhss
{
// For now, Client works as a "trusted key generator".
class Client
{
 public:
  using DecPlain = Poly[2];
  using SkPtr = std::unique_ptr<SecretKey>;
  using PkPtr = std::shared_ptr<Ciphertext>;

  Client(std::size_t hw = 64): sk_(nullptr), pk_(nullptr)
  { KeyGen(hw); }

  PkPtr GetPk() { return pk_; }

  /**
   * Performs LPR decryption on LPR ciphertext using BigPoly during scaling-and-rounding
   * @param[in] c LPR ciphertext to be decrypted
   * @param[out] m plaintext of single-precision magnitude
   */
  void Decrypt(const Ciphertext& c, Poly& m) const
  {
    // 1. c.b + c.a * s.a
    CRTPoly crt_dp;
    CRTPoly::MulMod(c.GetConstB(), sk_->GetConstB(), crt_dp);
    CRTPoly::AddModEqual(c.GetConstA(), crt_dp);
    BigNTT::ApplyInvNTT(crt_dp);
    // std::cout << crt_dp << std::endl;
    // 2. Scale and Rounding
    BigPoly bp;
    PolyConv::CRTToBigPoly(crt_dp, bp);
    // std::cout << bp << std::endl;
    
    /* Without symmetric rep, negative noise is considered to be a poly of super big magnitude */
    PolyConv::AsymmetricRepToSymmetricRepEqual(bp);
    // std::cout << "Centered Poly" << bp << std::endl;
    BigPoly::ScaleAndRound(bp, m);
  }

  void Decrypt(const Ciphertext& c, BigPoly& m) const
  {
    // 1. c.b + c.a * s.a
    CRTPoly crt_dp;
    CRTPoly::MulMod(c.GetConstB(), sk_->GetConstB(), crt_dp);
    CRTPoly::AddModEqual(c.GetConstA(), crt_dp);
    BigNTT::ApplyInvNTT(crt_dp);
    // std::cout << crt_dp << std::endl;
    // 2. Scale and Rounding
    PolyConv::CRTToBigPoly(crt_dp, m);
    /* Without symmetric rep, negative noise is considered to be a poly of super big magnitude */
    PolyConv::AsymmetricRepToSymmetricRepEqual(m);
    // std::cout << "Centered Poly" << bp << std::endl;
    BigPoly::ScaleAndRoundEqual(m);
  }

  /**
   * Performs secret sharing on the secret key
   * @param[in] original secret key (1, s) in NTT form
   * @param[out] sk0 a random element in (R_Q)^2
   * @param[out] sk1 sk - sk0
   */
  void ShareSk(SecretKey& sk0, SecretKey& sk1)
  {
    CRTSampler::SampleUniformModQ(sk0.GetB());
    CRTSampler::SampleUniformModQ(sk0.GetA());
    SecretKey::SubMod(*sk_, sk0, sk1);
  }

#if 0
  void ShareSk(SkPtr sk0, SkPtr sk1)
  {
    if (sk0 == nullptr)
       throw std::invalid_argument("sk0 is not initialized.");

    if (sk1 == nullptr)
       throw std::invalid_argument("sk1 is not initialized.");

    CRTSampler::SampleUniformModQ(sk0->GetB());
    CRTSampler::SampleUniformModQ(sk0->GetA());
    SecretKey::SubMod(*sk_, *sk0, *sk1);
  }
#endif
  /**
   * Performs secret sharing on the secret key
   * @param[in] original secret key (1, s) in NTT form
   * @param[out] sk0 (0, s0) \in (R_q)^2 for random s0 in NTT form
   * @param[out] sk1 (1, s - s0) \in (R_q)^2 
   */
  void ShareSkOpt(SecretKey& sk0, SecretKey& sk1)
  {
    CRTPoly::SetConstantNTT(1, sk0.GetB());
    CRTSampler::SampleUniformModQ(sk0.GetA());
    CRTPoly::SubMod(sk_->GetA(), sk0.GetA(), sk1.GetA());
  }

 private:
  /**
   * Performs a secret and public key pair under Ring-LWE assumption
   * @param[in] hm hamming weight of the secret key polynomial
   * @param[out] secret key and public key
   */
  void KeyGen(std::size_t hw)
  {
    sk_ = std::make_unique<SecretKey>();
    // 1. Create a secret key (1, s) in NTT domain
    CRTPoly::SetConstantNTT(1, sk_->GetB());
    CRTSampler::SampleUniformTernaryWithHW(hw, sk_->GetA());
    BigNTT::ApplyNTT(sk_->GetA());
    
    // 2. Create an encryption of zero under the secret key and publish it as a public key in NTT Domain
    // pk := (e - a * s , a)
    pk_ = std::make_shared<Ciphertext>();
    EncryptZero(*sk_, *pk_);
  }
  SkPtr sk_;
  PkPtr pk_;
};
} // namespace lhss