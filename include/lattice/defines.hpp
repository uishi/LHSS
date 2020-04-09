#pragma once

using u128 = unsigned __int128;
using i128 = __int128;

using UIntType = std::uint64_t;
using SIntType = std::int64_t;
using DoubleUIntType = u128;
using DoubleSIntType = i128;
using FloatType = double;
constexpr std::size_t kWordSize = 64;
constexpr std::size_t kWordSizeMinusOne = kWordSize - 1;

#include <gmpxx.h>