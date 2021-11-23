# Numerical methods 4 / Eigenvalues

Contains implementations of following methods for eigenvalue computation:

* QR-algorithm
* Modifiend reverse iteration method utilizing Rayleigh quotient
* Optional usage of shifts (selectable through QRMode)
* Optional usage of Heisenberg form (selectable by applying heisenberg_form())

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: MSVC v142
* Requires C++17 support

## Usage

Place config file of the following format into the same folder as executable:

* Line 1: INPUT_FILEPATH [value without whitespaces]
* Line 2: OUTPUT_FOLDER [value without whitespaces]
* Line 2: PRECISION

Refer to 'config.txt' as an example. Upon execution no furter inputs are required.

## Version history

* 00.04
    * Implemeted reverse iteration metohd
    * Implemented modifiend reverse iteration method utilizing Rayleigh quotient

* 00.03
    * Implemented computation of Hessenberg matrix
    * Added tests for QR-algorithm with Hessenberg form

* 00.02
    * Implemented shifts for QR-algorithm as an optional mode

* 00.01
    * Basic boilerplate
    * Implemented QR-algorithm

## License

This project is licensed under the MIT License - see the LICENSE.md file for details