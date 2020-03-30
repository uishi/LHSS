#pragma once

#include <cmath>
#include <random>

#include "params.hpp"
#include "poly.hpp"

namespace lhss
{
std::normal_distribution<double> d(0.0, 3.2);
// for ternary uniform sampling (pk gen and ctxt gen)
std::uniform_int_distribution<UIntType> z3_uniform(0, 2);
// for uniform sampling w/ hamming weight (sk gen)
std::uniform_int_distribution<UIntType> z2_uniform(0, 1);
struct SmallSampler
{
  using SmallPoly = Poly;

  // static std::random_device rd;

//   // @param poly  
//   static void SamplePolyFromUniformTernary(SmallPoly& p) {
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_int_distribution<> dis(-1, 1);
//     for (size_t i = 0; i != Params::n; ++i) {
//       int val = dis(gen);
//       if (val == -1) {
//         p[i] = Poly::Zq::modulus - 1;
//       } else {
//         p[i] = val;
//       }
//     }
//   }

//  @param poly  
  static void SamplePolyFromUniformMod(SmallPoly& p, UIntType modulus)
  {
    std::random_device rd("/dev/urandom");
    std::uniform_int_distribution<UIntType> dis(0, modulus);
    for (size_t i = 0; i != Params::n; ++i)
      p.coeffs[i] = dis(rd);
  }
};

struct CRTSampler
{
  // Samples gaussiann variables over R, then reduce it into modulo R_qi for i = 1,... ,k
  static void SamplePolyFromNormalDist(CRTPoly& p)
  {
    std::random_device rd("/dev/urandom");
    std::vector<UIntType> noise_poly(Params::n);
    std::vector<UIntType> signs(Params::n);
    std::size_t sample_cnt = 0;

    while (sample_cnt < Params::n)
    {
      std::int64_t r = d(rd);
      if (static_cast<UIntType>(std::abs(r)) > Params::Berr)
        continue;

      noise_poly[sample_cnt] = static_cast<UIntType>(r);
      signs[sample_cnt] = static_cast<UIntType>(r >> kWordSizeMinusOne);
      ++sample_cnt;
    }

    for (std::size_t i = 0; i != Params::num_moduli; ++i)
    {
      UIntType mod = Params::z_qi[i].modulus;
      for (std::size_t j = 0; j < Params::n; ++j)
      {
        p.small_polys[i].coeffs[j] = noise_poly[j];
         // if the noise is negative, adds the modulus to it
        p.small_polys[i].coeffs[j] += mod & signs[j];
      }
    }
  }

  // Samples gaussiann variables over R, followed by multiplying P to it.
  // Then reduce it into modulo R_qi for i = 1,... ,k
  static void SamplePolyFromNormalDistAndMultByConst(const std::vector<UIntType>& scales, CRTPoly& p)
  {
    std::random_device rd("/dev/urandom");
    for (std::size_t j = 0; j != Params::n; ++j)
    {
      std::int64_t r = d(rd);
      bool is_big = false;
      if (static_cast<UIntType>(std::abs(r)) > Params::Berr)
        is_big = true;
    
      for (std::size_t i = 0; i != Params::num_moduli; ++i)
      {
        // negative
        if (std::signbit(r))
        {
          if (is_big)
            p.small_polys[i].coeffs[j] = Params::z_qi[i].modulus - Params::Berr;
          else
            p.small_polys[i].coeffs[j] = Params::z_qi[i].modulus + static_cast<UIntType>(r);
        }
        else
        {
          if (is_big)
            p.small_polys[i].coeffs[j] = Params::Berr;
          else
            p.small_polys[i].coeffs[j] = static_cast<UIntType>(r);
        }
        Params::z_qi[i].MulModEqual(scales[i], p.small_polys[i].coeffs[j]);
      }
    }
  }

  static void SampleUniformTernaryWithHW(const std::size_t hw, CRTPoly& p)
  {
    std::random_device rd("/dev/urandom");
    std::uniform_int_distribution<std::size_t> zn_uniform(0, Params::n - 1);

    std::size_t non_zero_count = 0;
    std::size_t non_zero_coeff_id;

    while(non_zero_count != hw)
    {
      non_zero_coeff_id = zn_uniform(rd);
      if (p.small_polys[0].coeffs[non_zero_coeff_id] != 0)
        continue;
    
      UIntType val = z2_uniform(rd);
      for (std::size_t i = 0; i < Params::num_moduli; ++i)
      {
        if (val == 0) // 1
          p.small_polys[i].coeffs[non_zero_coeff_id] = 1;
        else  // -1
          p.small_polys[i].coeffs[non_zero_coeff_id] = Params::z_qi[i].modulus - 1;
      }
      ++non_zero_count;
    }
  }

  static void SampleUniformTernary(CRTPoly& p)
  {
    std::random_device rd("/dev/urandom");
    for (std::size_t j = 0; j != Params::n; ++j)
    {
      UIntType val = z3_uniform(rd);
      for (std::size_t i = 0; i != Params::num_moduli; ++i)
      {
        if (val == 1)  // 1
          p.small_polys[i].coeffs[j] = 1;
        else if (val == 2) // -1
          p.small_polys[i].coeffs[j] = Params::z_qi[i].modulus - 1;
      }
    }
  }

  static void SampleUniformModQ(CRTPoly& p)
  {
    for (std::size_t i = 0; i != Params::num_moduli; ++i)
      SmallSampler::SamplePolyFromUniformMod(p.small_polys[i], Params::z_qi[i].modulus);
  }
};
} // namespace lhss