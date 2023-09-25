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
#include "../LSQ.h"

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
    cout << "Test para el metodo de minimos cuadrados" << endl;
    cout << "Test de la matriz transpuesta" << endl;
    cout << endl;

    std::vector<double> simpleData = {1.0, 3.0, -1.0, 13.0, 4.0, -1.0, 1.0, 9.0, 2.0, 4.0, 3.0, -6.0};
    Matrix<double> testMatrix(3, 4, simpleData);

    cout << "Matriz Original" << endl;
    testMatrix.printMatrix();
    cout << endl;

    Matrix<double> testMatrixT = testMatrix.Transpose();
    cout << "La matriz transpuesta es: " << endl;
    testMatrixT.printMatrix();

    cout << endl << "Test de minimos cuadrados" << endl;
    std::vector<double> Xdata = {1.0, 1.0, 1.0, 2.0, 1.0, 3.0};
    Matrix<double> X(3, 2, &Xdata);

    std::vector<double> Ydata = {2.0, 4.0, 4.0};
    Vector<double> y(Ydata);

    cout << "Matriz X = " << endl;
    X.printMatrix();
    cout << endl;
    cout << "Y vector y = " << endl;
    PrintVector(y);
    cout << endl;

    Vector<double> result;
    int test = LSQ(X, y, result);

    cout << "El mejor beta^ = " << endl;
    PrintVector(result);
    cout << endl;
}