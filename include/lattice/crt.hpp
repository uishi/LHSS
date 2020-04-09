#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include "defines.hpp"
#include "util.hpp"

namespace lhss
{
struct CRT
{
  std::vector<UIntType> qis;
  std::size_t nummod;
  std::vector<UIntType> q_div_qi_invs_mod_qis;
  mpz_class q;// for debugging
  std::vector<mpz_class> q_div_qis; // for debugging

  CRT(){}
  
  // Offline
  CRT(const std::vector<UIntType>& moduli)
   : qis(moduli), nummod(qis.size())
  {
    if (nummod <= 1)
      throw std::invalid_argument("# CRT moduli must be > 1");
  
    q_div_qi_invs_mod_qis.reserve(nummod);
    PreComputeCRTInvs();

    q_div_qis.reserve(nummod);
    PreComputeQandQdivQis();
  }

  void Reconstruct(const std::vector<UIntType>& rems, mpz_class& val)
  {
    if (rems.size() != qis.size())
      throw std::invalid_argument("# rems must be equal to # CRT moduli");
  
    val = 0;

    mpz_class rem = rems[0];
    rem *= q_div_qi_invs_mod_qis[0];
    rem *= q_div_qis[0];
    val += rem;
    for (std::size_t i = 1; i < nummod; ++i)
    {
      rem = rems[i] * q_div_qi_invs_mod_qis[i];
      rem *= q_div_qis[i];
      val += rem;
    }  
    val %= q;
  }

  void ReconstructWithIntermediateModReduce(const std::vector<UIntType>& rems, mpz_class& val)
  {
    if (rems.size() != qis.size())
      throw std::invalid_argument("# rems must be equal to # CRT moduli");

    val = 0;

    mpz_class rem = rems[0];
    rem *= q_div_qi_invs_mod_qis[0];
    rem *= q_div_qis[0];
    val += rem;
    for (std::size_t i = 1; i < nummod; ++i) 
    {
      rem = rems[i];
      rem *= q_div_qi_invs_mod_qis[i];
      rem %= qis[i];
      rem *= q_div_qis[i];
      val += rem;
    }  
    val %= q;
  }

#if 0
  void ComputeReconstructRem(const std::vector<UIntType>& rems, mpz_class& val)
  {
    if (rems.size() != qis.size())
      throw std::invalid_argument("# rems must be equal to # CRT moduli");
  
    val = 0;

    mpz_class rem = rems[0];
    rem *= q_div_qi_invs_mod_qis[0];
    rem *= q_div_qis[0];
    val += rem;
    for (std::size_t i = 1; i < nummod; ++i)
    {
      rem = rems[i] * q_div_qi_invs_mod_qis[i];
      rem *= q_div_qis[i];
      val += rem;
    }
    val /= q;
  }
#endif

  void Set(const std::vector<UIntType>& moduli)
  {
    qis = moduli;
    nummod = qis.size();
    if (nummod <= 1)
      throw std::invalid_argument("# CRT moduli must be > 1");

    q_div_qi_invs_mod_qis.reserve(nummod);
    PreComputeCRTInvs();

    q_div_qis.reserve(nummod);
    PreComputeQandQdivQis();
  }

 private:
  // Valid for many moduli (k > 2)
  void PreComputeCRTInvs()
  {
    std::vector<UIntType> qi_inv_mod_qj(nummod*nummod, 0);
    for (std::size_t i = 0; i < nummod; ++i) 
    {
      for (std::size_t j = 0; j < nummod; ++j) 
      {
        if (i == j) continue;
        MultInverse(qis[j], qis[i], qi_inv_mod_qj[i*nummod + j]);
      }
      
      std::uint64_t q_div_qi_inv_mod_qi;
      if (i == 0)
      {
        q_div_qi_inv_mod_qi = qi_inv_mod_qj[1];
        // FIXME
        for (std::size_t j = 2; j < nummod; ++j)
          MultModEqualNaive(qi_inv_mod_qj[j], qis[0], q_div_qi_inv_mod_qi);
      }
      else
      { 
        q_div_qi_inv_mod_qi = qi_inv_mod_qj[i * nummod];
        for (std::size_t j = 1; j < nummod; ++j)
        {
          if (i == j)
            continue;

          MultModEqualNaive(qi_inv_mod_qj[i * nummod + j], qis[i], q_div_qi_inv_mod_qi);
        }
      }
      q_div_qi_invs_mod_qis.push_back(q_div_qi_inv_mod_qi);
    }
  }

  void PreComputeQandQdivQis()
  {
    q = 1;
    for (std::size_t i = 0; i < nummod; ++i)
      q *= qis[i];

    for (std::size_t i = 0; i < nummod; ++i)
      q_div_qis.push_back(q/qis[i]);
  }
};
}