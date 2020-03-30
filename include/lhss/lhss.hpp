#pragma once

#include <stdexcept>

#include "lattice.hpp"
#include "rlwe_ops.hpp"
#include "rlwe_ops_sym.hpp"
#include "secretkey.hpp"
#include "evaluator.hpp"
#include "client.hpp"

namespace lhss
{
using SkShare = SecretKey;
using SkPtr = std::unique_ptr<SecretKey>;
using PkPtr = std::shared_ptr<Ciphertext>;

void LHSSInit(
 const std::vector<UIntType>& plain_mods,
 const std::vector<UIntType>& ctxt_mods,
 const std::size_t deg
)
{
  Params::Init(plain_mods, ctxt_mods, deg);
}

// Failure probability per multiplication is <= 2^-40.
void Setup(const std::size_t max_bound_message_bit)
{
  /**  log(P) = 62, log(Q) = 144, N = 2^13, security >= 128
   * - Size of P: {62}
   * - Size of Q: {62, 41, 41}
   */
  if (max_bound_message_bit == 1)
    LHSSInit({4611686018427322369}, {4611686018427322369, 2199023190017, 2199022927873}, 1 << 13);
  /**  log(P) = 77, log(Q) = 174, N = 2^13, security >= 80
   * - P = {q1, q2}; size: {38, 39}
   * - Q = {q1, q2, q3, q4}; size: {38, 39, 48, 49}
   */
  else if (max_bound_message_bit == 16)
    LHSSInit({274877562881, 549755731969}, {274877562881, 549755731969, 281474976694273, 562949952847873}, 1 << 13);
  /**  log(P) = 93, log(Q) = 206, N = 2^13, security >= 80
   * - P = {q1, q2};         size: {46, 47}
   * - Q = {q1, q2, q3, q4}; size: {46, 47, 56, 57}
   */
  else if (max_bound_message_bit == 32)
    LHSSInit({70368743669761, 140737488273409}, {70368743669761, 140737488273409, 72057594037616641, 144115188075593729}, 1 << 13);
  /**  log(P) = 125, log(Q) = 270, N = 2^13, security >= 80
   * - P = {q1, q2, q3};             size: {41, 42, 42}
   * - Q = {q1, q2, q3, q4, q5, q6}; size: {41, 42, 42, 48, 48, 49}
   */
  else if (max_bound_message_bit == 64)
    LHSSInit({2199023190017, 4398045708289, 4398046150657}, {2199023190017, 4398045708289, 4398046150657, 281474976694273, 281474976546817, 562949952847873}, 1 << 13);

  else
    throw std::invalid_argument("The maximum magnitude (in bits) must be either 1, 16, 32, or 64...");
}
} // namespace lhss