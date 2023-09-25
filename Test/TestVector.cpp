#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "../Vector.h"
#include "../Matrix.h"

using namespace std;

// A simple function to print a matrix to stdout.
template <class T>
void PrintVector(Vector<T> vector)
{
	int dims = vector.getNumDim();
	for (int i = 0; i < dims; ++i){
    cout << std::fixed << std::setprecision(3) << vector.getElement(i) << endl;
	}
    cout << endl;
}

int main(){
    cout << "Test para la clase Vector" << endl;
    cout << "Vector 1" << endl;
    std::vector<double> vectorData1;
    vectorData1.push_back(1.0);
    vectorData1.push_back(2.0);
    vectorData1.push_back(3.0);
    Vector<double> vector1(vectorData1);
    PrintVector(vector1);
    cout << "Vector 2" << endl;
    std::vector<double> vectorData2;
    vectorData2.push_back(3.0);
    vectorData2.push_back(2.0);
    vectorData2.push_back(1.0);
    Vector<double> vector2(vectorData2);
    PrintVector(vector2);
    cout << "Vector 1 - Vector 2 = " << endl;
    PrintVector(vector1 - vector2);
    cout << "Vector 1 + Vector 2 = " << endl;
    PrintVector(vector1 + vector2);
    cout << "Vector 1 * 2.0 = " << endl;
    PrintVector(vector1 * 2.0);
    cout << "2 * Vector 1 = " << endl;
    PrintVector(2.0 * vector1);
    cout << "Vector e(1)" << endl;
    std::vector<double> vectorData3;
    vectorData3.push_back(1.0);
    vectorData3.push_back(0.0);
    vectorData3.push_back(0.0);
    Vector<double> vector3(vectorData3);
    PrintVector(vector3);
    cout << "Vector e(2)" << endl;
    std::vector<double> vectorData4;
    vectorData4.push_back(0.0);
    vectorData4.push_back(1.0);
    vectorData4.push_back(0.0);
    Vector<double> vector4(vectorData4);
    PrintVector(vector4);
    cout << "Calculamos el producto vectoria de e(1) * e(2) = e(3)" << endl;
    PrintVector(Vector<double>::cross(vector3, vector4));
    cout << "Calculamos el producto escalar de Vector 1 * e(2)" << endl;
    cout << Vector<double>::dot(vector1, vector4) << endl;

    cout << "Calculamos el producto escalar de e(1) * e(2) = 0" << endl;
    cout << Vector<double>::dot(vector3, vector4) << endl;

    cout << "Norma del Vector 1: ";
    cout << vector1.norm() << endl;

    cout << "Calculo del vector normalizado de Vector 1 (funcion con retorno)" << endl;
    PrintVector(vector1.Normalized());
    cout << "Calculo del vector normalizado del Vector 2" << endl;
    vector2.Normalize();
    PrintVector(vector2);
    cout << "Multiplicacion de Matriz * Vector 1" << endl;
    double dataMatrix[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix<double> matrix(3, 3, dataMatrix);
    cout << "Matriz: " << endl;
    matrix.printMatrix();
    Vector<double> result = matrix * vector1;
    cout << "Resultado: " << endl;
    PrintVector(result);
}