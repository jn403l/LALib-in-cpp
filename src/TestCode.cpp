/* code to test the maths functions contained in the library*/

#include <iostream>
#include <string>
#include <math.h>

#include "qbMatrix.h"

using namespace std;

// a simple function to print a matrix to stdout
template<class T>
void PrintMatrix(qbMatrix2<T> matrix) {
  int nRows = matrix.GetNumRows();
  int nCols = matrix.GetNumCols();
  for (int row = 0; row < nRows; ++row) {
    for (int col = 0; col < nCols; ++col) {
      cout << matrix.GetElement(row, col) << " ";
    }
    cout << endl;
  }
}

int main() {
  cout << "Code to test qbMatrix2" << endl;

  // **************************************************
  // create an instance of the qbmatrix2 class
  // this will contain a simple 2D 3x3 matrix
  double simpleData[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  qbMatrix2<double> testMatrix(3, 4, simpleData);

  // extract and print the elements of test matrix
  cout << endl << "******************************" << endl;
  cout << "3x4 matrix test (testMatrix)." << endl;
  PrintMatrix(testMatrix);
  // **************************************************
	// Test element retrieval
	cout << endl << "********************" << endl;
	cout << "Element (0,0) = " << testMatrix.GetElement(0,0) << endl;
	cout << "Element (1,0) = " << testMatrix.GetElement(1,0) << endl;
	cout << "Element (2,0) = " << testMatrix.GetElement(2,0) << endl;
	cout << "Element (0,1) = " << testMatrix.GetElement(0,1) << endl;
	cout << "Element (1,1) = " << testMatrix.GetElement(2,1) << endl;
	cout << "Element (5,5) = " << testMatrix.GetElement(5,5) << endl;
  // **************************************************
	// Test matrix multiplication
	cout << endl << "********************" << endl;
	cout << "Test matrix multiplication." << endl;
	double simpleData2[12] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
	qbMatrix2<double> testMatrix2(4, 3, simpleData2);
	cout << "4x3 matrix (testMatrix2)" << endl;
	PrintMatrix(testMatrix2);
	cout << "Multiplication (testMatrix * testMatrix2) result:" << endl;
	qbMatrix2<double> multTest1 = testMatrix * testMatrix2;
	PrintMatrix(multTest1);

	cout << endl << "********************" << endl;
	cout << "testMatrix2 * testMatrix:" << endl;
	qbMatrix2<double> multTest2 = testMatrix2 * testMatrix;
	PrintMatrix(multTest2);

	cout << endl << "********************" << endl;
	cout << "Test multiplication of column vector by matrix." << endl;
	double columnData[3] = {1.5, 2.5, 3.5};
	double squareData[9] = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
	qbMatrix2<double> testColumn(3, 1, columnData);
	qbMatrix2<double> squareMatrix(3, 3, squareData);
	cout << "Column vector = " << endl;
	PrintMatrix(testColumn);
	cout << "Square matrix = " << endl;
	PrintMatrix(squareMatrix);
	cout << "Column vector * Square matrix = " << endl;
	PrintMatrix(testColumn * squareMatrix);
	cout << "Square matrix * Column vector = " << endl;
	PrintMatrix(squareMatrix * testColumn);
	cout << "Square matrix + 1.0 = " << endl;
	qbMatrix2<double> squareMatrix2 = squareMatrix + 1.0;
	PrintMatrix(squareMatrix2);
	cout << "(Square matrix + 1.0) + Column vector = " << endl;
	PrintMatrix(squareMatrix2 + testColumn);

  // **************************************************
	// Test equality operator
	cout << endl << "********************" << endl;
	cout << "Test equality operator." << endl;
	cout << "testMatrix == testMatrix2: " << (testMatrix == testMatrix2) << endl;
	cout << "testMatrix2 == testMatrix: " << (testMatrix2 == testMatrix) << endl;
  cout << "(Let testMatrix3 = testMatrix)" << endl;
  qbMatrix2<double> testMatrix3 = testMatrix;
	cout << "testMatrix == testMatrix3: " << (testMatrix == testMatrix3) << endl;
	cout << "testMatrix3 == testMatrix: " << (testMatrix3 == testMatrix) << endl;

  // **************************************************
	// Test matrix addition by scalar
	cout << endl << "********************" << endl;
	cout << "Test addition by scalar." << endl;
  cout << "testMatrix + 2.0 = " << endl;
  PrintMatrix(testMatrix+2.0);
  cout << endl;
  cout << "2.0 + testMatrix = " << endl;
  PrintMatrix(2.0+testMatrix);

  // **************************************************
	// Test matrix subtraction by scalar
	cout << endl << "********************" << endl;
	cout << "Test addition by scalar." << endl;
  cout << "testMatrix - 2.0 = " << endl;
  PrintMatrix(testMatrix-2.0);
  cout << endl;
  cout << "2.0 - testMatrix = " << endl;
  PrintMatrix(2.0-testMatrix);

  // **************************************************
	// Test matrix multiplication by scalar
	cout << endl << "********************" << endl;
	cout << "Test multiplication by scalar." << endl;
  cout << "Test * 2.0 = " << endl;
  PrintMatrix(testMatrix*2.0);
  cout << endl;
  cout << "2.0 * testMatrix = " << endl;
  PrintMatrix(2.0*testMatrix);

  return 0;
}
