#include <stdexcept>
#include <tuple>

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>

#include "tmatrix.hpp"
#include "qr_algorithm.hpp"
#include "hessenberg_form.hpp"
#include "reverse_iteration.hpp"



constexpr auto CONFIG_PATH = "config.txt";

void separator() {
	std::cout << "\n#########################";
}

// Parse config and return input/output filepaths
// @return 1 => input filepath
// @return 2 => output filepath
// @return 3 => precision
std::tuple<std::string, std::string, double> parse_config() {
	std::ifstream inConfig(CONFIG_PATH);

	if (!inConfig.is_open()) throw std::runtime_error("Could not open config file");

	std::string outputFolder;
	std::string inputPath;
	double precision;

	std::string dummy; // used to skip comments

	inConfig
		>> dummy >> inputPath
		>> dummy >> outputFolder
		>> dummy >> precision;

	return { inputPath, outputFolder, precision };
}

DMatrix parse_matrix(const std::string &filepath) {
	std::ifstream inFile(filepath);

	if (!inFile.is_open())
		throw std::runtime_error("ERROR: Could not open file <" + filepath + ">");

	// Read size, construct matrix
	size_t n;
	inFile >> n;
	DMatrix A(n, n);

	// Read matrix elements in 1 direct-indexation loop
	for (size_t i = 0; i < A.total_size(); ++i) inFile >> A(i);

	return A;
}

DMatrix run_qr_algorithm(const DMatrix &A, double precision, QRMode mode, QRNote note, const std::string &outputPath) {
	separator();
	std::cout
		<< "\n>>> Computing eigenvalues..."
		<< "\n>>> Method    -> QR algorithm"
		<< "\n>>> Precision -> " << precision
		<< "\n>>> Shifts    -> " << (static_cast<bool>(mode) ? "ON" : "OFF");

	auto [eigenvals, iterations, operations] = QR_algorithm(A, precision, mode, note);

	std::cout
		<< "\n>>> Iterations = " << iterations
		<< "\n>>> Operations = " << operations << "\n";
	eigenvals.print();

	// Save result into the file
	std::ofstream outFile(outputPath + '[' + (static_cast<bool>(mode) ? "ON" : "OFF") + ']' + ".txt");
	
	outFile
		<< "Iterations: " << iterations
		<< "\nOperations: " << operations
		<< "\nEigenvalues:";

	for (size_t i = 0; i < eigenvals.total_size(); ++i) outFile << '\n' << eigenvals(i);

	return eigenvals;
	// File is closed automatically
}

void run_reverse_iteration(const DMatrix &A, const DMatrix &eigenvals, double precision, const std::string &outputPath) {
	separator();
	std::cout
		<< "\n>>> Computing eigenvectors..."
		<< "\n>>> Method    -> Reverse iteration"
		<< "\n>>> Precision -> " << precision;

	// Save result into the file 1 by 1
	std::ofstream outFile(outputPath + "[eigenvecs].txt");

	for (size_t k = 0; k < eigenvals.total_size(); ++k) {
		const auto eigenval = eigenvals(k);
		const auto [eigenvec, iterations, residual] = reverse_iteration(A, eigenval, precision);

		std::cout
			<< "\n>>> Eigenvector [" << k << "]:"
			<< '\n' << TMATRIX_PRINT_INDENT << "lambda_" << k << "   = " << eigenval
			<< '\n' << TMATRIX_PRINT_INDENT << "iterations = " << iterations
			<< '\n' << TMATRIX_PRINT_INDENT << "residual   = " << residual;

		outFile
			<< "\n[" << k << "]:"
			<< "\nlambda_" << k << "   -> " << eigenval
			<< "\niterations -> " << iterations
			<< "\nresidual   -> " << residual
			<< "\neigenvector:";
		for (size_t i = 0; i < eigenvec.total_size(); ++i) outFile << '\n' << eigenvec(i);
		outFile << '\n';
	}

	std::cout << "\n>>>Done! Result saved to " << outputPath << "[eigenvecs].txt";
}


void run_rayleigh_reverse_iteration(const DMatrix &A, double precision, const std::string &outputPath) {
	separator();
	std::cout
		<< "\n>>> Computing eigenvalues and eigenvectors..."
		<< "\n>>> Method    -> Rayleigh reverse iteration"
		<< "\n>>> Precision -> " << precision;

	const auto N = A.rows();

	// Save result into the file 1 by 1
	std::ofstream outFile(outputPath + "[rayleigh].txt");

	DMatrix initialEstimate(N, 1);

	for (size_t k = 0; k < N; ++k) {
		// Choose initial estimate to your liking
		fill(initialEstimate, 0.);
		initialEstimate(k) = 1.;

		std::cout << "\n>>> Initial estimate [" << k << "]:\n";
		initialEstimate.print();

		const auto [value, vector, iterations, residual] = rayleigh_reverse_iteration(A, initialEstimate, precision);

		std::cout
			<< '\n' << TMATRIX_PRINT_INDENT << "lambda_" << k << "   = " << value
			<< '\n' << TMATRIX_PRINT_INDENT << "iterations = " << iterations
			<< '\n' << TMATRIX_PRINT_INDENT << "residual   = " << residual
			<< '\n' << TMATRIX_PRINT_INDENT << "vector_" << k << "   =\n";

		vector.print();

		outFile
			<< "\n[" << k << "]:"
			<< "\n estimate   ->\n";
		for (size_t i = 0; i < initialEstimate.total_size(); ++i) outFile << '\n' << initialEstimate(i);

		outFile
			<< "\nlambda_" << k << "    -> " << value
			<< "\niterations  -> " << iterations
			<< "\nresidual    -> " << residual
			<< "\neigenvector ->\n";
		for (size_t i = 0; i < vector.total_size(); ++i) outFile << '\n' << vector(i);
		outFile << '\n';
	}

	std::cout << "\n>>>Done! Result saved to " << outputPath << "[rayleigh].txt\n";
}


int main(int* argc, char** argv) {
	try {
		// Parse config
		std::cout << ">>> Parsing config...\n";
		const auto [inputPath, outputPath, precision] = parse_config();

		std::cout << ">>> Input path = " << inputPath << "\n>>> Output folder = " << outputPath << "\n";

		// Parse matrix
		std::cout << ">>> Parsing matrix...\n";
		const auto A = parse_matrix(inputPath);

		separator();
		std::cout << "\n>>> A(" << A.rows() << ", " << A.cols() << "):\n";
		A.print();

		// QR algorithm
		const auto eigenvals = run_qr_algorithm(A, precision, QRMode::SHIFTS_ON, QRNote::REGULAR, outputPath + "[Default]");
		run_qr_algorithm(A, precision, QRMode::SHIFTS_OFF, QRNote::REGULAR, outputPath + "[Default]");

		// Get Hesenberg matrix
		const auto hessenbergForm = hessenberg_form(A);

		separator();
		std::cout << "\n>>> HessenbergForm(" << hessenbergForm.rows() << ", " << hessenbergForm.cols() << "):\n";
		hessenbergForm.print();

		// QR algorithm with Hessenberg matrix
		run_qr_algorithm(hessenbergForm, precision, QRMode::SHIFTS_ON, QRNote::HESSENBERG, outputPath + "[Hessenberg]");
		run_qr_algorithm(hessenbergForm, precision, QRMode::SHIFTS_OFF, QRNote::HESSENBERG, outputPath + "[Hessenberg]");

		// Find eigenvectors through reverse iteration
		run_reverse_iteration(A, eigenvals, precision, outputPath);

		// Find eigenvalues and eigenvectors through
		run_rayleigh_reverse_iteration(A, precision, outputPath);
	}
	// If caught any errors, show error message
	catch (const std::runtime_error& err) {
		std::cerr << "\nRUNTIME EXCEPTION -> " << err.what() << std::endl;
	}
	catch (...) {
		std::cerr << "\nCAUGHT UNKNOWN EXCEPTION" << std::endl;
	}

	return 0;
}