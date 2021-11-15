# Numerical methods 4 / Eigenvalues

Contains implementations of following methods for eigenvalue computation:

* QR-algorithm

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: MSVC v142
* Requires C++17 support

## Usage

Place config file of the following format into the same folder as executable:

* Line 1: INPUT_FILEPATH [value without whitespaces]
* Line 2: OUTPUT_FOLDER [value without whitespaces]

Refer to 'config.txt' as an example. Upon execution no furter inputs are required.

## Version history

* 00.01
    * Basic boilerplate
    * Implemented QR-algorithm

## License

This project is licensed under the MIT License - see the LICENSE.md file for details