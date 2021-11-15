#include <stdexcept>
#include <tuple>

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>

#include "tmatrix.hpp"
#include "qr_algorithm.hpp"
#include "hessenberg_form.hpp"



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

void run_qr_algorithm(const DMatrix &A, double precision, QRMode mode, const std::string &outputFolder) {
	separator();
	std::cout
		<< "\n>>> Method    -> QR algorithm"
		<< "\n>>> Precision -> " << precision
		<< "\n>>> Shifts    -> " << (static_cast<bool>(mode) ? "ON" : "OFF");

	auto [eigenvals, iterations] = QR_algorithm(A, precision, mode);

	std::cout << "\n>>> Eigenvalues computed in " << iterations << " iterations:\n";
	eigenvals.print();

	/// SAVE TO THE FILE
}


int main(int* argc, char** argv) {
	try {
		// Parse config
		std::cout << ">>> Parsing config...\n";
		const auto [inputPath, outputFolder, precision] = parse_config();

		std::cout << ">>> Input path = " << inputPath << "\n>>> Output folder = " << outputFolder << "\n";

		// Parse matrix
		std::cout << ">>> Parsing matrix...\n";
		const auto A = parse_matrix(inputPath);

		separator();
		std::cout << "\n>>> A(" << A.rows() << ", " << A.cols() << "):\n";
		A.print();

		// QR algorithm
		run_qr_algorithm(A, precision, QRMode::SHIFTS_ON, outputFolder);
		run_qr_algorithm(A, precision, QRMode::SHIFTS_OFF, outputFolder);

		// Get Hesenberg matrix
		const auto hessenbergForm = hessenberg_form(A);

		separator();
		std::cout << "\n>>> HessenbergForm(" << hessenbergForm.rows() << ", " << hessenbergForm.cols() << "):\n";
		hessenbergForm.print();

		// QR algorithm with Hessenberg matrix
		run_qr_algorithm(hessenbergForm, precision, QRMode::SHIFTS_ON, outputFolder);
		run_qr_algorithm(hessenbergForm, precision, QRMode::SHIFTS_OFF, outputFolder);
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