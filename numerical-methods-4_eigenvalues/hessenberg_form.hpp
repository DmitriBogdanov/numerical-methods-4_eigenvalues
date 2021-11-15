#pragma once

#include "tmatrix.hpp"



void primitive_rotation_matrix(DMatrix &T, const DMatrix &A, size_t i, size_t j) {
	const double inverseDenominator = 1. / std::sqrt(sqr(A[i][i - 1]) + sqr(A[j][i - 1]));
	const double alpha = A[i][i - 1] * inverseDenominator;
	const double beta = A[j][i - 1] * inverseDenominator;

	// Build T_ij
	fill(T, 0.);
	fill_diagonal(T, 1.);
	if (i >= j) {
		T[i][j] = -beta;
		T[j][i] = beta;
	}
	else {
		T[i][j] = beta;
		T[j][i] = -beta;
	}
	
	T[i][i] = alpha;
	T[j][j] = alpha;
}

inline DMatrix hessenberg_form(const DMatrix &A) {
	const auto N = A.rows();

	auto HES = A;

	DMatrix T(N, N);
	DMatrix inverseT(N, N);
	DMatrix temp(N, N);

	for (size_t i = 2; i < N; ++i)
		for (size_t j = 0; j < i - 1; ++j) {
			// Build T_ij
			primitive_rotation_matrix(T, HES, j + 1, i);

			// Get T_ij^-1
			inverseT = T;
			transpose_square(inverseT);

			// HES = T_ij HES T_ij^-1
			multiply(temp, T, HES);
			multiply(HES, temp, inverseT);
		}

	return HES;
}