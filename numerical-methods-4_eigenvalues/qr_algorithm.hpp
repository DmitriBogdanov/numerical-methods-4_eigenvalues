#pragma once

#include <tuple>

#include "tmatrix.hpp"



enum class QRMode {
	SHIFTS_ON = 1,
	SHIFTS_OFF = 0
};

enum class QRNote {
	REGULAR,
	HESSENBERG
};

// Shift all diagonal elements by value
void diagonal_shift(DMatrix &matrix, double shiftValue) {
	for (size_t i = 0; i < matrix.rows(); ++i) matrix[i][i] += shiftValue;
}


// @return 1 => column of eigenvals
// @return 2 => number of iterations
// @return 2 => number of operations
inline std::tuple<DMatrix, uint, uint> QR_algorithm(const DMatrix &matrix, double precision, QRMode mode, QRNote note) {
	const auto N = matrix.rows();

	uint iterations = 0;
	uint operations = (note == QRNote::HESSENBERG) ? static_cast<uint>(4. * cube(N) - 1.5 * sqr(N) - 2.5 * N) : 0;

	DMatrix Ak(matrix), Anext(matrix);

	DMatrix Qk(N, N), Rk(N, N);

	double totalShift = 0.;
	size_t m = 0;

	while (m < N) {
		// Shift diagonal by sigma = A[m][m]
		// - only if correct mode is enabled
		if (mode == QRMode::SHIFTS_ON) {
			const double sigma = Ak[m][m];
			diagonal_shift(Ak, -sigma);
			totalShift += sigma;
		}

		// A_k+1 = R_k Q_k + sigma E
		QR_decompose(Qk, Rk, Ak);
		multiply(Ak, Rk, Qk);

		operations += (note == QRNote::HESSENBERG) ? sqr(N - m) : cube(N - m) * 8 / 3;
			// QR decomposition takes 8/3 n^3 operations
			// But for hessenberg only n^2
		operations += cube(N - m); // Multiplying matrices
		
		// Check if values to the left are small enough so we can go to the next 'm'
		double maxValue = 0.;
		for (size_t j = 0; j < m; ++j) maxValue = std::max(maxValue, std::abs(Ak[m][j]));

		if (maxValue < precision) ++m;

		++iterations;
	}

	// Store eigenvals as a column, account for shifts
	DMatrix eigenvals(N, 1);
	for (size_t i = 0; i < N; ++i) eigenvals(i) = Ak[i][i] + totalShift;

	return { eigenvals, iterations, operations };
}