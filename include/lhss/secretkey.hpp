#pragma once

#include "poly.hpp"

namespace lhss
{
struct SecretKey
{
  CRTPoly elements[2];

  const CRTPoly& GetConstB() const { return elements[0]; }

  const CRTPoly& GetConstA() const { return elements[1]; }

  CRTPoly& GetB() { return elements[0]; }

  CRTPoly& GetA() { return elements[1]; }

  static void AddMod(const SecretKey& s0, const SecretKey& s1, SecretKey& s_add)
  {
    CRTPoly::AddMod(s0.GetConstB(), s1.GetConstB(), s_add.GetB());
    CRTPoly::AddMod(s0.GetConstA(), s1.GetConstA(), s_add.GetA());
  }

  static void AddMod(const SecretKey& s, SecretKey& s_add)
  {
    CRTPoly::AddMod(s.GetConstB(), s_add.GetB());
    CRTPoly::AddMod(s.GetConstA(), s_add.GetA());
  }


  static void SubMod(const SecretKey& sk0, const SecretKey& sk1, SecretKey& sk_sub)
  {
    CRTPoly::SubMod(sk0.GetConstB(), sk1.GetConstB(), sk_sub.GetB());
    CRTPoly::SubMod(sk0.GetConstA(), sk1.GetConstA(), sk_sub.GetA());
  }

  static void ApplyAurotmophismNTT(const SecretKey& s_in, const std::size_t amt, SecretKey& s_out)
  {
    CRTPoly::AutomorphismNTT(s_in.GetConstB(), amt, s_out.GetB());
    CRTPoly::AutomorphismNTT(s_in.GetConstA(), amt, s_out.GetA());
  }
};
} // namespace lhss