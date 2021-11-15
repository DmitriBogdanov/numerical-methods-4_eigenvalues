#include "tmatrix.hpp"



// @return 1 => column of eigenvals
inline DMatrix QR_algorith(const DMatrix &matrix, double precision) {
	const auto N = matrix.rows();

	DMatrix Ak(matrix), Anext(matrix);

	DMatrix Qk(N, N), Rk(N, N);

	uint iter = 0;

	while (!is_close_to_upper_triangular(Ak, precision)) {
		++iter;

		QR_decompose(Qk, Rk, Ak);

		// A_k+1 = R_k Q_k
		multiply(Ak, Rk, Qk);
	}

	// Store eigenvals as a column
	DMatrix eigenvals(N, 1);
	for (size_t i = 0; i < N; ++i) eigenvals(i) = Ak[i][i];

	return eigenvals;
}