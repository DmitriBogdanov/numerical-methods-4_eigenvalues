#pragma once

#include "tmatrix.hpp"



// @return 1 => eigenvector
// @return 2 => number of iterations
// @return 2 => residual
inline std::tuple<DMatrix, uint, double> reverse_iteration(const DMatrix &A, double lambda, double precision) {
	const auto N = A.rows();

	// Denote C = A - lambda E
	DMatrix C = A;
	for (size_t i = 0; i < N; ++i) C[i][i] -= lambda;

	/// Prepare linear solver
	///DMatrix Q(N, N), R(N, N);
	///QR_decompose(Q, R, C);

	// Residuals
	DMatrix residualVector(N, 1);
	double residual = INF;

	// Initial estimate is (1, 0, ... , 0)
	DMatrix X0(N, 1);
	fill(X0, 0.);
	X0(0) = 1.;

	DMatrix X(N, 1);

	uint iterations = 0;

	while (true) {
		// Find X_k+1 by solving linear system C*X_k+1 = X_k
		X = gauss_solve(C, X0);
		///QR_solve(X, Q, R, X0);
		

		// Normalize solution, rotate it if necessary
		const auto factor = 1. / eucledian_norm(X) * (X(0) > 0 ? 1. : -1.);
		for (size_t i = 0; i < N; ++i) X(i) *= factor;

		// Find residual
		multiply(residualVector, C, X);
		residual = eucledian_norm(residualVector);

		// Stopping criteria
		++iterations;
		if (vector_difference_norm(X, X0) < precision) break;

		// X0 now X
		X0 = X;
	}

	return { X, iterations, residual };
}



// @return 1 => eigenvalue
// @return 2 => eigenvector
// @return 3 => iterations
// @return 4 => residual
inline std::tuple<double, DMatrix, uint, double> rayleigh_reverse_iteration(const DMatrix &A, const DMatrix &initialEstimate, double precision) {
	const auto N = A.rows();

	DMatrix temp(N, 1);
	DMatrix C(N, N);

	double lambda;

	// Residuals
	DMatrix residualVector(N, 1);
	double residual = INF;
	
	// Normalize X0 just in case
	auto X0 = initialEstimate;
	const auto factor = 1. / eucledian_norm(X0);
	for (size_t i = 0; i < N; ++i) X0(i) *= factor;

	DMatrix X(N, 1);

	uint iterations = 0;

	while (true) {
		// 1) Get lambda_k
		multiply(temp, A, X0);
		lambda = dot_product(temp, X0); // respects ||X0|| = 1

		// Set C = A - lambda_k E
		C = A;
		for (size_t i = 0; i < N; ++i) C[i][i] -= lambda;

		// 2) Get X_k+1
		X = gauss_solve(C, X0);

		// 3) Normalize X_k+1
		const auto factor = 1. / eucledian_norm(X) * (X(0) > 0 ? 1. : -1.);
		for (size_t i = 0; i < N; ++i) X(i) *= factor;

		// Find residual
		multiply(residualVector, C, X);
		residual = eucledian_norm(residualVector);

		// Stopping criteria
		++iterations;
		if (vector_difference_norm(X, X0) < precision) break;

		// X now X0
		X0 = X;
	}

	return { lambda, X, iterations, residual };
}