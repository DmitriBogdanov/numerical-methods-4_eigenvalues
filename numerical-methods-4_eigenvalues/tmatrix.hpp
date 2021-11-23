#pragma once

#include <vector> // internal storage uses std::vector
#include <algorithm> // std::max
#include <type_traits> // std::is_arithmetic<T> is used to ensure numeric underlying type
#include <iostream> // console output
#include <iomanip> // console output format

#include "math_helpers.hpp"



// Config
constexpr size_t TMATRIX_MAX_PRINT_SIZE = 8;
constexpr std::streamsize TMATRIX_PRINT_ALIGNMENT = 14;
constexpr const char *TMATRIX_PRINT_INDENT = "   ";



// # TMatrix #
// Lightweight vector wrapper for neat matr
// - Contiguous data storage
// - standart indexation -> A[i][j]
// - for row/column matrices 1D indexation can be used -> A(k)
template<typename T>
class TMatrix {
	// Template compiles only if T is numeric
	static_assert(std::is_arithmetic<T>::value, "T must be numeric");

	size_t _rows;
	size_t _cols;
	std::vector<T> _data;

public:
	TMatrix() :
		_rows(0),
		_cols(0)
	{}

	TMatrix(size_t rows, size_t cols) :
		_rows(rows),
		_cols(cols),
		_data(rows * cols)
	{}

	// 1D indexation (k)
	T& operator() (size_t k) { return _data[k]; }
	const T& operator() (size_t k) const { return _data[k]; }

	// 2D indexation [i][j]
	T* operator[] (size_t i) { return _data.data() + i * _cols; } // allows [i][j] style of indexation
	const T* operator[] (size_t i) const { return _data.data() + i * _cols; }

	// Getters
	size_t rows() const { return _rows; };
	size_t cols() const { return _cols; };
	size_t total_size() const { return _rows * _cols; }

	void print() const {
		// If any of matrix dimensions is larger than that console output is supressed
		if (_rows > TMATRIX_MAX_PRINT_SIZE || _cols > TMATRIX_MAX_PRINT_SIZE) {
			std::cout << TMATRIX_PRINT_INDENT << "[ matrix output supressed due to large size ]\n";
			return;
		}

		for (size_t i = 0; i < _rows; ++i) {
			std::cout << TMATRIX_PRINT_INDENT << "[ ";

			const auto IxC = i * _cols;
			for (size_t j = 0; j < _cols; ++j) std::cout << std::setw(TMATRIX_PRINT_ALIGNMENT) << _data[IxC + j] << " ";

			std::cout << " ]\n";
		}
	}

	// Cubic norm
	T norm() const {
		// Cubic norm => max_i { SUM_j |matrix[i][j]|}
		T maxSum(0);

		// Go over each row calculating SUM_j |matrix[i][j]|, select max_i
		for (size_t i = 0; i < _rows; ++i) {
			T sum(0);

			const auto IxC = i * _cols;
			for (size_t j = 0; j < _cols; ++j) sum += std::abs(_data[IxC + j]);

			maxSum = std::max(maxSum, sum);
		}

		return maxSum;
	}

	/// Octahedral norm (uncomment if necessary)
	//T norm_octahedral() const {
	//	// Octahedral norm => max_i { SUM_j |matrix[i][j]|}
	//	T maxSum(0);

	//	// Go over each row calculating SUM_j |matrix[i][j]|, select max_i
	//	for (size_t j = 0; j < _cols; ++j) {
	//		T sum(0);

	//		for (size_t i = 0; i < _rows; ++i) sum += std::abs(_data[i * _cols + j]);

	//		maxSum = std::max(maxSum, sum);
	//	}

	//	return maxSum;
	//}
};



// # Template shortcuts #
using DMatrix = TMatrix<double>;
using FMatrix = TMatrix<float>;
using IMatrix = TMatrix<int>;



// # Matrix operations #
// - To discourage unnecessary allocation, all matrix operations are implemented as functions 
//   that store their result into an externally allocated matrix passed by reference &dest
template<typename T>
void fill(TMatrix<T> &dest, T value) {
	for (size_t k = 0; k < dest.total_size(); ++k) dest(k) = value;
}

template<typename T>
void fill_diagonal(TMatrix<T> &dest, T value) {
	for (size_t k = 0; k < dest.rows(); ++k) dest[k][k] = value;
}

// matrix + matrix
template<typename T>
void add(TMatrix<T> &dest, const TMatrix<T> &src1, const TMatrix<T> &src2) {
	for (size_t k = 0; k < dest.total_size(); ++k) dest(k) = src1(k) + src2(k);
}

// matrix - matrix
template<typename T>
void substract(TMatrix<T> &dest, const TMatrix<T> &src1, const TMatrix<T> &src2) {
	for (size_t k = 0; k < dest.total_size(); ++k) dest(k) = src1(k) - src2(k);
}

// matrix * matrix
template<typename T>
void multiply(TMatrix<T> &dest, const TMatrix<T> &src1, const TMatrix<T> &src2) {
	fill(dest, static_cast<T>(0));

	for (size_t i = 0; i < src1.rows(); ++i)
		for (size_t k = 0; k < src1.cols(); ++k)
			for (size_t j = 0; j < src2.cols(); ++j)
				dest[i][j] += src1[i][k] * src2[k][j];
				// note that naive loop order would be [i]->[j]->[k], swapping [k] and [j]
				// loops reduces the number of cache misses since we access contiguously
				// stored elements in the inner-most loop
}

// matrix * scalar
template<typename T>
void multiply(TMatrix<T> &dest, const TMatrix<T> &src, T value) {
	for (size_t k = 0; k < dest.total_size(); ++k) dest(k) = src(k) * value;
}

// vector * vector
template<typename T>
T dot_product(TMatrix<T> &vec1, const TMatrix<T> &vec2) {
	T sum(0);

	for (size_t k = 0; k < vec1.total_size(); ++k) sum += vec1(k) * vec2(k);

	return sum;
}

// Returns eucledian vector norm
template<typename T>
T eucledian_norm(const TMatrix<T> &vector) {
	T sum(0);

	for (size_t i = 0; i < vector.total_size(); ++i) sum += sqr(vector(i));

	return std::sqrt(sum);
}

// Returns ||A-B|| where A, B are column-vectors
template<typename T>
T vector_difference_norm(const TMatrix<T> &src1, const TMatrix<T> &src2) {
	T norm(0);
	for (size_t i = 0; i < src1.total_size(); ++i) norm = std::max(norm, std::abs(src1(i) - src2(i)));

	return norm;
}

template<typename T>
bool is_close_to_upper_triangular(const TMatrix<T> &matrix, T precision) {
	T maxFound(0);

	for (size_t i = 0; i < matrix.rows(); ++i)
		for (size_t j = 0; j < i; ++j)
			maxFound = std::max(maxFound, std::abs(matrix[i][j]));

	return maxFound < precision;
}

// Solves tridiagonal system and stores solution into the &dest
// - saves on copying but modifies A, b during computation
template<typename T>
TMatrix<T>& tridiagonal_solve_in_place(TMatrix<T> &A, TMatrix<T> &b) {
	const auto N = A.rows();

	// Forward elimination
	for (size_t i = 0; i < N - 1; ++i) {
		// Eliminate element below [i][i]
		const T factor = A[i + 1][0] / A[i][1];

		A[i + 1][0] -= A[i][1] * factor;
		A[i + 1][1] -= A[i][2] * factor;
		b(i + 1) -= b(i) * factor;
	}

	// Backward elimination
	for (size_t i = N - 1; i > 0; --i) {
		// Eliminate element above [i][i]
		const T factor = A[i - 1][2] / A[i][1];

		A[i - 1][2] -= A[i][1] * factor;
		b(i - 1) -= b(i) * factor;
	}

	// Normalize diagonal
	for (size_t i = 0; i < N; ++i) {
		const T factor = 1. / A[i][1];

		A[i][1] *= factor;
		b(i) *= factor;
	}

	// Now 'b' holds the solution
	return b;
}

// Transposes square A
template<typename T>
void transpose_square(TMatrix<T> &A) {
	const auto N = A.rows();

	for (size_t i = 0; i < N; ++i)
		for (size_t j = i; j < N; ++j)
			std::swap(A[i][j], A[j][i]);
}

template<typename T>
DMatrix gauss_solve(const TMatrix<T> &inpA, const TMatrix<T> &inpB) {
	const auto N = inpA.rows();

	// Concate matrices A = inpA|inpB
	TMatrix<T> A(N, N + 1);

	for (size_t i = 0; i < N; ++i)
		for (size_t j = 0; j < N; ++j)
			A[i][j] = inpA[i][j];

	for (size_t i = 0; i < N; ++i)
		A[i][N] = inpB[i][0];

	// Forward elimination
	for (size_t k = 0; k < A.rows(); ++k) {
		// Select leading row and swap it with current row if needed
		auto leadingRow = k;
		for (size_t i = k; i < A.rows(); ++i)
			if (std::abs(A[i][k]) > std::abs(A[leadingRow][k]))
				leadingRow = i;

		// Swap row with leading row
		if (leadingRow != k)
			for (size_t j = 0; j < A.cols(); ++j)
				std::swap(A[leadingRow][j], A[k][j]);

		// If the diagonal element of a leading row is zero, matrix is singular => throw
		if (is_zero(A[k][k]))
			throw std::runtime_error("ERROR: Could not solve the system with singular matrix");

		// Set diagonal element to 1 by multiplying its row
		const T factor = static_cast<T>(1) / A[k][k];
		for (size_t j = k; j < A.cols(); ++j) A[k][j] *= factor;

		// Substract current row from the following ones completing the step
		for (size_t i = k + 1; i < A.rows(); ++i) {
			const T firstElement = -A[i][k];

			for (size_t j = k; j < A.cols(); ++j) A[i][j] += A[k][j] * firstElement;
		}
	}

	// Backward elimination
	// Go backwards until our matrix is diagonal
	for (size_t k = A.rows() - 1; k > 0; --k)
		for (size_t i = 1; i <= k; ++i) {
			const T factor = -A[k - i][k];
			for (size_t j = k; j <= A.cols(); ++j) A[k - i][j] += A[k][j] * factor;
		}

	// Return solution as a column
	TMatrix<T> solution(N, 1);
	for (size_t i = 0; i < A.rows(); ++i) solution[i][0] = A[i][A.cols() - 1];

	return solution;
}

// Returns QR decomposition of A
// - outputs into Q, R
template<typename T>
void QR_decompose(TMatrix<T> &Q, TMatrix<T> &R, const TMatrix<T> &A) {
	const auto N = A.rows();

	// Set Q to identity matrix
	fill(Q, 0.);

	for (size_t i = 0; i < N; ++i)
		Q[i][i] = 1;

	// Copy original matrix into R
	R = A;

	// Givens rotations go brrr
	T c(0);
	T s(0);
	for (size_t i = 0; i < N - 1; ++i)
		for (size_t j = i + 1; j < N; ++j)
			if (!is_zero(R[j][i])) {
				const T inverseNorm = static_cast<T>(1) / std::sqrt(sqr(R[i][i]) + sqr(R[j][i]));
				c = R[i][i] * inverseNorm;
				s = R[j][i] * inverseNorm;

				for (size_t j_ = 0; j_ < N; ++j_) {
					const T tempR = c * R[i][j_] + s * R[j][j_];
					R[j][j_] = -s * R[i][j_] + c * R[j][j_];
					R[i][j_] = tempR;

					const T tempQ = c * Q[i][j_] + s * Q[j][j_];
					Q[j][j_] = -s * Q[i][j_] + c * Q[j][j_];
					Q[i][j_] = tempQ;
				}
			}

	transpose_square(Q);
}

// Solve linear system using existing QR-decomposition
// - stores solution in 'dest'
template<typename T>
void QR_solve(TMatrix<T> &dest, const TMatrix<T> &inpQ, const TMatrix<T> &inpR, const TMatrix<T> &b) {
	const auto N = inpQ.rows();

	// Get Q^T
	auto R = inpR;
	auto QT = inpQ;
	transpose_square(QT);

	// Get right side f = Q^T * b
	multiply(dest, QT, b);

	// Solve through Gauss backwards elimination
	for (size_t k = 0; k < N; ++k) {
		const T factor = 1. / R[k][k];
		for (size_t j = k; j < N; ++j) R[k][j] *= factor;
		dest[k][0] *= factor;
	}

	for (size_t k = N - 1; k > 0; --k) {
		for (size_t i = 1; i <= k; ++i) {
			const T factor = R[k - i][k];
			R[k - i][k] -= R[k][k] * factor;
			dest[k - i][0] -= dest[k][0] * factor;
		}
	}
}