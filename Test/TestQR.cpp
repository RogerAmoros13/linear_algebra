
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <sstream>
#include <vector>
#include <random>
#include <fstream>

#include "../Matrix.h"
#include "../Vector.h"
#include "../QR.h"

using namespace std;

template <class T>
void PrintVector(Vector<T> inputVector)
{
	int nRows = inputVector.getNumDim();
	for (int row = 0; row<nRows; ++row)
  {
  cout << std::fixed << std::setprecision(3) << inputVector.getElement(row) << endl;
	}
}

int main(){
	cout << "************************************************" << endl;
	cout << "Test de la descomposicion QR" << endl;
	cout << "************************************************" << endl;
	cout << endl;
	{
		cout << "Test con matriz 3 x 3" << endl;

		std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
		Matrix<double> testMatrix (3, 3, simpleData);

		testMatrix.printMatrix();

		cout << endl;
		cout << "Calculo descomposicion QR..." << endl;

		Matrix<double> Q (3, 3);
		Matrix<double> R (3, 3);

		int status = QR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;

		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;

		cout << "Q x R = " << endl;
		(Q*R).printMatrix();
	}

	{
		cout << "Test con matriz 4 x 4" << endl;

		std::vector<double> simpleData = {1.0, 5.0, 3.0, 4.0, 7.0, 8.0, 2.0, 9.0, 7.0, 3.0, 2.0, 1.0, 9.0, 3.0, 5.0, 7.0};
		Matrix<double> testMatrix(4, 4, simpleData);

		testMatrix.printMatrix();

		cout << endl;
		cout << "Calculo descomposicion QR..." << endl;

		Matrix<double> Q (4, 4);
		Matrix<double> R (4, 4);

		int status = QR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;

		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;

		cout << "Q x R = " << endl;
		(Q*R).printMatrix();
	}

	{
		cout << "Test con matriz 5 x 5" << endl;

		std::vector<double> simpleData = {2, 6, 4, 6, 8, 6, 7, 9, 7, 9, 2, 3, 6, 3, 5, 6, 1, 1, 5, 5, 3, 5, 6, 5, 6};
		Matrix<double> testMatrix(5, 5, simpleData);

		testMatrix.printMatrix();

		cout << endl;
		cout << "Calculo descomposicion QR..." << endl;

		Matrix<double> Q (4, 4);
		Matrix<double> R (4, 4);

		int status = QR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;

		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;

		cout << "Q x R = " << endl;
		(Q*R).printMatrix();
	}

	{
		cout << "Test con matriz 5 x 5 float" << endl;

		std::vector<double> simpleData = {
			8.662634278267483, 2.3440981169711796, 3.414158790068152, 9.819959485632891, 9.812414578216162,
			4.8096369839436495, 7.743133259609277, 9.871217856632036, 7.100783013043249, 8.127838524397976,
			1.3468248609110365, 1.3120774834063536, 9.607366488550678, 2.852679282078192, 8.087038227451359,
			7.556075051454403, 5.80117852857823, 3.550189544341768, 3.7807047754393994, 7.934423413357392,
			2.866445996919499, 7.125441061546031, 4.53141730712106, 4.297092147605687, 2.5126585000174146
		};
		Matrix<double> testMatrix(5, 5, simpleData);

		testMatrix.printMatrix();

		cout << endl;
		cout << "Calculo descomposicion QR..." << endl;

		Matrix<double> Q (4, 4);
		Matrix<double> R (4, 4);

		int status = QR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;

		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;

		cout << "Q x R = " << endl;
		(Q*R).printMatrix();
	}
}
