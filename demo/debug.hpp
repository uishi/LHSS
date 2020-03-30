// should be included after `lhss.hpp`
#pragma once

using namespace lhss;
using SkShare = SecretKey;

static void DebugShare(const SkShare& s0, const SkShare& s1, bool is_all = false) {
  SkShare t0 = s0;
  SkShare t1 = s1;
  SecretKey debug_m;
  SecretKey::AddMod(t0, t1, debug_m);

  BigNTT::ApplyInvNTT(debug_m.GetB());

  BigPoly debug_message[2];
  PolyConv::CRTToBigPoly(debug_m.GetConstB(), debug_message[0]);
  std::cout << "\t  Output  (x):";
  BigPoly::ModReduceEqual(Params::P, debug_message[0]);
  PolyConv::AsymmetricRepToSymmetricRepEqual(debug_message[0]);
  BigPoly::ShowSomeCoeff(debug_message[0], 0, 10, std::cout);
  std::cout << std::endl;

  if (is_all)
  {
    BigNTT::ApplyInvNTT(debug_m.GetA());
    PolyConv::CRTToBigPoly(debug_m.GetConstA(), debug_message[1]);
    std::cout << "\t  Output(x*s):";
    BigPoly::ModReduceEqual(Params::P, debug_message[1]);
    PolyConv::AsymmetricRepToSymmetricRepEqual(Params::P, debug_message[1]);
    BigPoly::ShowSomeCoeff(debug_message[1], 0, 10, std::cout);
    std::cout << std::endl;
  }
}

static void DebugCtxt(const SecretKey& sk, const HSSCtxt& ct, bool is_all = false) {
  BigPoly m_org[2];
  DecryptOrgSk(sk, ct.enc_m, m_org[0]);
  PolyConv::AsymmetricRepToSymmetricRepEqual(Params::P, m_org[0]);
  std::cout << "\t  Dec(Enc(m)):";
  BigPoly::ShowSomeCoeff(m_org[0], 0, 10, std::cout);
  std::cout << std::endl;

  if (is_all)
  {
    DecryptOrgSk(sk, ct.enc_m_times_s, m_org[1]);
    PolyConv::AsymmetricRepToSymmetricRepEqual(Params::P, m_org[1]);
    std::cout << "\t  Dec(Enc(m*s)):";
    BigPoly::ShowSomeCoeff(m_org[1], 0, 10, std::cout);
    std::cout << std::endl;
  }
}