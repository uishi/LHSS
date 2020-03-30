#pragma once

#include "poly.hpp"

namespace lhss
{
struct Ciphertext
{
  CRTPoly elements[2];

  CRTPoly& GetB() { return elements[0]; }
  CRTPoly& GetA() { return elements[1]; }
  const CRTPoly& GetConstB() const { return elements[0]; }
  const CRTPoly& GetConstA() const { return elements[1]; }

  static void AddMod(const Ciphertext& c1, const Ciphertext& c2, Ciphertext& c_add)
  {
    CRTPoly::AddMod(c1.GetConstB(), c2.GetConstB(), c_add.GetB());
    CRTPoly::AddMod(c1.GetConstA(), c2.GetConstA(), c_add.GetA());
  }

  static void SubMod(const Ciphertext& c1, const Ciphertext& c2, Ciphertext& c_sub)
  {
    CRTPoly::SubMod(c1.GetConstB(), c2.GetConstB(), c_sub.GetB());
    CRTPoly::SubMod(c1.GetConstA(), c2.GetConstA(), c_sub.GetA());
  }
};

struct HSSCtxt
{  
  Ciphertext enc_m;
  Ciphertext enc_m_times_s;

  static void AddMod(const HSSCtxt& c1, const HSSCtxt& c2, HSSCtxt& c_add)
  { 
    Ciphertext::AddMod(c1.enc_m, c2.enc_m, c_add.enc_m);
    Ciphertext::AddMod(c1.enc_m_times_s, c2.enc_m_times_s, c_add.enc_m_times_s);
  }

  static void Negate(HSSCtxt& c)
  {
    CRTPoly::NegateEqual(c.enc_m.GetB());
    CRTPoly::NegateEqual(c.enc_m.GetA());
    CRTPoly::NegateEqual(c.enc_m_times_s.GetB());
    CRTPoly::NegateEqual(c.enc_m_times_s.GetA());
  }
};

struct HSSCiphertext
{  
  CRTPoly elements[4];

  CRTPoly& GetB0() { return elements[0]; }
  CRTPoly& GetA0() { return elements[1]; }
  CRTPoly& GetB1() { return elements[2]; }
  CRTPoly& GetA1() { return elements[3]; }
  const CRTPoly& GetConstB0() const { return elements[0]; }
  const CRTPoly& GetConstA0() const { return elements[1]; }
  const CRTPoly& GetConstB1() const { return elements[2]; }
  const CRTPoly& GetConstA1() const { return elements[3]; }

  static void AddMod(const HSSCiphertext& c1, const HSSCiphertext& c2, HSSCiphertext& c_add)
  { 
    CRTPoly::AddMod(c1.GetConstB0(), c2.GetConstB0(), c_add.GetB0());
    CRTPoly::AddMod(c1.GetConstA0(), c2.GetConstA0(), c_add.GetA0());
    CRTPoly::AddMod(c1.GetConstB1(), c2.GetConstB1(), c_add.GetB1());
    CRTPoly::AddMod(c1.GetConstA1(), c2.GetConstA1(), c_add.GetA1());
  }

  static void AddConstantNTTModEqual(UIntType constant, HSSCiphertext& c)
  {
    // TODO FIXME do this ahead of time
    CRTPoly scaled_constant;
    Poly m;
    m.SetDenseSameElem(constant);
    CRTPoly::MulSmallPolyAndCRTConstant(m, Params::q_div_p_mod_qi, scaled_constant);

    CRTPoly::AddModEqual(scaled_constant, c.GetB0());
    CRTPoly::AddModEqual(scaled_constant, c.GetA1());
  }

  static void SubMod(const HSSCiphertext& c1, const HSSCiphertext& c2, HSSCiphertext& c_sub)
  {  
    CRTPoly::SubMod(c1.GetConstB0(), c2.GetConstB0(), c_sub.GetB0());
    CRTPoly::SubMod(c1.GetConstA0(), c2.GetConstA0(), c_sub.GetA0());
    CRTPoly::SubMod(c1.GetConstB1(), c2.GetConstB1(), c_sub.GetB1());
    CRTPoly::SubMod(c1.GetConstA1(), c2.GetConstA1(), c_sub.GetA1());
  }

  static void NegateEqual(HSSCiphertext& c)
  {
    CRTPoly::NegateEqual(c.GetB0());
    CRTPoly::NegateEqual(c.GetA0());
    CRTPoly::NegateEqual(c.GetB1());
    CRTPoly::NegateEqual(c.GetA1());
  }
};


struct LWECiphertext
{
  CRTInt b;   // Z_Q
  CRTPoly a;  // (Z_Q)^n
};

#if 0
struct CipherConv
{
  static void ConvRLWEToLWE(Ciphertext& c_src, LWECiphertext& c_dest)
  {
  }

  // // supports once Ciphertext data structur is changed.
  // static void ConvRLWEToLWEInplace(Ciphertext& c_src)
  // {    
  // }
};
#endif

} // namespace lhss