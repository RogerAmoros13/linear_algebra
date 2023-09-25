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
#include "../EIG.h"

using namespace std;

template <class T>
void PrintVector(Vector<T> inputVector)
{
	int nRows = inputVector.getNumDim();
	for (int row = 0; row<nRows; ++row)
  {
  cout << std::fixed << std::setprecision(6) << inputVector.getElement(row) << endl;
	}
}

int main(){
  cout << "***************************************" << endl;
  cout << "Test de valores y vectores propios" << endl;
  cout << "Power Iteration Method" << endl;
  cout << "***************************************" << endl;
  {
    cout << endl;
    cout << "Matriz 2x2: " << endl;
    std::vector<double> matrixData = {1.0, 2.0, 3.0, 4.0};
    Matrix<double> matrix(2, 2, matrixData);
    matrix.printMatrix();

    cout << endl;
    cout << "Calculamos el vector y el valor propio" << endl;
    double eigenValue;
    Vector<double> eigenVector;
    int test = EIG_PIt(matrix, eigenValue, eigenVector);

    cout << "Vector propio: " << endl;
    PrintVector(eigenVector);
    cout << endl << "Valor propio: " << eigenValue << endl;

    cout << endl;
    cout << "Matriz 3x3: " << endl;
    std::vector<double> matrixData2 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix<double> matrix2(3, 3, matrixData2);
    matrix2.printMatrix();

    cout << endl;
    cout << "Calculamos el vector y el valor propio" << endl;
    double eigenValue2;
    Vector<double> eigenVector2;
    int test2 = EIG_PIt(matrix2, eigenValue2, eigenVector2);

    cout << "Vector propio: " << endl;
    PrintVector(eigenVector2);
    cout << endl << "Valor propio: " << eigenValue2 << endl;
  }

  cout << "***************************************" << endl;
  cout << "Test de valores y vectores propios" << endl;
  cout << "Inverse-Power Iteration Method." << endl;
  cout << "***************************************" << endl;
  {
    cout << "Testing with a simple 3x3 matrix:" << endl;
    std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
    Matrix<double> testMatrix(3, 3, simpleData);
    testMatrix.printMatrix();
    cout << endl;

    std::vector<double> eigenValues = {1.5962551, -0.37253087, 0.02627577};
    Vector<double> eigenVector(3);

    for (auto currentValue : eigenValues){
      cout << "Vector propio estimado del valor propio " << currentValue << " = " << endl;
      int returnStatus = InvPIt<double>(testMatrix, currentValue, eigenVector);
      PrintVector(eigenVector);

      if (returnStatus == EIG_MAXITERATIONSEXCEEDED){
        cout << "*** Maximo de iteraciones alcanzado ***" << endl;
      }
      cout << endl;
    }
  }

  cout << "****************************************" << endl;
  cout << "Test de valores y vectores propios" << endl;
  cout << "Funcion IsSymmetric en la clase de Matrix" << endl;
  cout << "*****************************************" << endl;
  {
    cout << "Test con matriz simetrica" << endl;
    std::vector<double> simpleData = {4.0, -7.0, 6.0, -7.0, -2.0, 13.0, 6.0, 13.0, 5.0};
    Matrix<double> testMatrix(3, 3, simpleData);
    testMatrix.printMatrix();
    cout << endl;

    cout << "La matriz es simetrica: ";
    if (testMatrix.isSymmetric()){
      cout << "Verdadero!" << endl;
    }else{
      cout << "Falso :(" << endl;
    }
  }

  {
    cout << "Test con matriz no-simetrica" << endl;
    std::vector<double> simpleData = {5.0, -4.0, 6.0, -7.0, -2.0, 13.0, 6.0, 13.0, 5.0};
    Matrix<double> testMatrix(3, 3, simpleData);
    testMatrix.printMatrix();
    cout << endl;

    cout << "La matriz es simetrica: ";
    if (testMatrix.isSymmetric()){
      cout << "Verdadero!" << endl;
    }else{
      cout << "Falso :(" << endl;
    }
  }

  cout << endl << endl;
  cout << "****************************************" << endl;
  cout << "Test de valores y vectores propios" << endl;
  cout << "Funcion IsSymmetric en la clase de Matrix" << endl;
  cout << "*****************************************" << endl;

  {
    cout << "Test con matriz 3x3 no-simetrica: " << endl;
    std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
    Matrix<double> testMatrix (3, 3, simpleData);
    testMatrix.printMatrix();
    cout << endl;

    std::vector<double> eigenValues;
    int returnStatus = EigQR(testMatrix, eigenValues);

    if (returnStatus == EIG_MAXITERATIONSEXCEEDED){
      cout << "*** Maximo de iteraciones excedido ***" << endl;
    }

    if (returnStatus == EIG_MATRIXNOTSYMMETRIC){
      cout << "*** Matriz no simetrica ***" << endl;
    }

    cout << "Los valores propios estimados son: " << endl;
    for (auto currentValue : eigenValues){
      cout << setprecision(3) << currentValue << " ";
    }
    cout << endl << endl;
  }

  {
    cout << "Test con matriz 3x3 simetrica: " << endl;
    std::vector<double> simpleData = {6.0, 5.5, -1.0, 5.5, 1.0, -2.0, -1.0, -2.0, -3.0};
    Matrix<double> testMatrix (3, 3, simpleData);
    testMatrix.printMatrix();
    cout << endl;

    std::vector<double> eigenValues;
    int returnStatus = EigQR(testMatrix, eigenValues);

    if (returnStatus == EIG_MAXITERATIONSEXCEEDED){
      cout << "*** Maximo de iteraciones excedido ***" << endl;
    }

    if (returnStatus == EIG_MATRIXNOTSYMMETRIC){
      cout << "*** Matriz no simetrica ***" << endl;
    }

    Vector<double> eigenVector(3);

    cout << "Los valores propios estimados son: " << endl;
    for (auto currentValue : eigenValues){
      cout << "Vector propio estimado para el valor propio " << currentValue << " = " << endl;
      int returnStatus = InvPIt<double>(testMatrix, currentValue, eigenVector);
      PrintVector(eigenVector);
      if (returnStatus == EIG_MAXITERATIONSEXCEEDED){
        cout << "*** Maximo de iteraciones excedido ***" << endl;
      }
      cout << endl;
    }
    cout << endl << endl;
  }

  {
    cout << "Test con matriz que deberia tener valores propios complejos" << endl;
    std::vector<double> simpleData = {4.0, -6.0, 8.0, 7.0, 9.0, -5.0, 9.0, -6.0, -4.0};
    Matrix<double> testMatrix(3, 3, simpleData);
    testMatrix.printMatrix();
    cout << endl;

    std::vector<double> eigenValues;
    int returnStatus = EigQR(testMatrix, eigenValues);

    if (returnStatus == EIG_MAXITERATIONSEXCEEDED){
      cout << "*** Maximo de iteraciones excedido ***" << endl;
    }

    if (returnStatus == EIG_MATRIXNOTSYMMETRIC){
      cout << "*** Matriz no simetrica ***" << endl;
    }

    cout << "Los valores propios estimados son: " << endl;
    for (auto currentValue : eigenValues){
      cout << std::setprecision(6) << currentValue << " ";
    }

    cout << endl << endl;

    Vector<double> eigenVector(3);

    // cout << "Los valores propios estimados son: " << endl;
    // for (auto currentValue : eigenValues){
    //   cout << "Vector propio estimado para el valor propio " << currentValue << " = " << endl;
    //   int returnStatus = InvPIt<double>(testMatrix, currentValue, eigenVector);
    //   PrintVector(eigenVector);
    //   if (returnStatus == EIG_MAXITERATIONSEXCEEDED){
    //     cout << "*** Maximo de iteraciones excedido ***" << endl;
    //   }
    //   cout << endl;
    // }
    // cout << endl << endl;
  }
}