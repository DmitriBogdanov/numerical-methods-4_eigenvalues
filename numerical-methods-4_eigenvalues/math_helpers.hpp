#pragma once

#include <algorithm>


using uint = unsigned int;
using sint = int;

constexpr double INF = std::numeric_limits<double>::infinity();

template<typename T>
constexpr bool is_zero(T value) {
	constexpr T EPSILON = static_cast<T>(1e-16);

	return std::abs(value) < EPSILON;
}

// -1, 0, 1 signum
template<typename T>
constexpr T sign(T value) {
	return (T(0) < value) - (value < T(0));
}

template<typename T>
constexpr T sqr(T value) {
	return value * value;
}

template<typename T>
constexpr T cube(T value) {
	return value * value * value;
}
