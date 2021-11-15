#include <stdexcept>
#include <tuple>

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>

#include "tmatrix.hpp"
#include "qr_algorithm.hpp"



constexpr auto CONFIG_PATH = "config.txt";

// Parse config and return input/output filepaths
// @return 1 => input filepath
// @return 2 => output filepath
std::tuple<std::string, std::string> parse_config() {
	std::ifstream inConfig(CONFIG_PATH);

	if (!inConfig.is_open()) throw std::runtime_error("Could not open config file");

	std::string outputFolder;
	std::string inputPath;

	std::string dummy; // used to skip comments

	inConfig
		>> dummy >> inputPath
		>> dummy >> outputFolder;

	return { inputPath, outputFolder };
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



int main(int* argc, char** argv) {
	try {
		// Parse config
		std::cout << ">>> Parsing config...\n";
		const auto [inputPath, outputFolder] = parse_config();

		std::cout << ">>> Input path = " << inputPath << "\n>>> Output folder = " << outputFolder << "\n";

		// Parse matrix
		std::cout << ">>> Parsing matrix...\n";
		const auto A = parse_matrix(inputPath);

		std::cout << ">>> A(" << A.rows() << ", " << A.cols() << "):\n";
		A.print();

		QR_algorith(A, 1e-6);

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