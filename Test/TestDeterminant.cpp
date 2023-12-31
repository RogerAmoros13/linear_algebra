#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "../Matrix.h"

using namespace std;

// A simple function to print a matrix to stdout.
template <class T>
void PrintMatrix(Matrix<T> matrix)
{
	int nRows = matrix.getNumRows();
	int nCols = matrix.getNumCols();
	for (int row = 0; row<nRows; ++row)
  {
	  for (int col = 0; col<nCols; ++col)
    {
	    cout << std::fixed << std::setprecision(3) << matrix.GetElement(row, col) << "  ";
    }
	cout << endl;
	}
}

int main()
{

	cout << "Testing implementation of determinant calculation." << endl;
	cout << endl;

	cout << "Generate a test matrix." << endl;
	double testData[9] = {2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};
	Matrix<double> testMatrix(3, 3, testData);
	PrintMatrix(testMatrix);
	cout << endl;

	cout << "Extract sub-matrix for element (0,0)" << endl;
	Matrix<double> minor1 = testMatrix.FindSubMatrix(0,0);
	PrintMatrix(minor1);
	cout << endl;

	cout << "Extract sub-matrix for element (0,1)" << endl;
	Matrix<double> minor2 = testMatrix.FindSubMatrix(0,1);
	PrintMatrix(minor2);
	cout << endl;

	cout << "Extract sub-matrix for element (0,2)" << endl;
	Matrix<double> minor3 = testMatrix.FindSubMatrix(0,2);
	PrintMatrix(minor3);
	cout << endl;

	cout << "Extract sub-matrix for element (1,1)" << endl;
	Matrix<double> minor4 = testMatrix.FindSubMatrix(1,1);
	PrintMatrix(minor4);
	cout << endl;

	cout << "Test with a larger matrix." << endl;
	double testData2[25] =
		{2.0, 3.0, 4.0, 5.0, 6.0,
		 1.0, 2.0, 3.0, 4.0, 5.0,
		 9.0, 5.0, 3.0, 2.0, 6.0,
		 2.0, 4.0, 6.0, 5.0, 1.0,
		 1.0, 7.0, 5.0, 2.0, 3.0};
	Matrix<double> testMatrix2(5, 5, testData2);
	PrintMatrix(testMatrix2);
	cout << endl;

	cout << "Extract sub-matrix for element (0,0)" << endl;
	Matrix<double> minor5 = testMatrix2.FindSubMatrix(0,0);
	PrintMatrix(minor5);
	cout << endl;

	cout << "Extract sub-matrix for element (0,1)" << endl;
	Matrix<double> minor6 = testMatrix2.FindSubMatrix(0,1);
	PrintMatrix(minor6);
	cout << endl;

	cout << "Extract sub-matrix for element (0,2)" << endl;
	Matrix<double> minor7 = testMatrix2.FindSubMatrix(0,2);
	PrintMatrix(minor7);
	cout << endl;

	cout << "Extract sub-matrix for element (1,1)" << endl;
	Matrix<double> minor8 = testMatrix2.FindSubMatrix(1,1);
	PrintMatrix(minor8);
	cout << endl;

	cout << "Test determinant of 3x3 matrix:" << endl;
	cout << testMatrix.Determinant() << endl;
	cout << endl;

	cout << "Test determinant of 5x5 matrix:" << endl;
	cout << testMatrix2.Determinant() << endl;
	cout << endl;

	cout << "Test determinant of a singular matrix:" << endl;
	double testData3[9] =
		{1.0, 1.0, 1.0,
		 0.0, 1.0, 0.0,
		 1.0, 0.0, 1.0};
	Matrix<double> testMatrix3(3, 3, testData3);
	PrintMatrix(testMatrix3);
	cout << endl;
	cout << "Determinant = " << testMatrix3.Determinant() << endl;
	cout << endl;

	return 0;
}