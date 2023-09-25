#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <random>

#include "../Matrix.h"
#include "../Vector.h"
#include "../LinSolve.h"

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

// A simple function to print a vector to stdout.
template <class T>
void PrintVector(Vector<T> inputVector)
{
	int nRows = inputVector.getNumDim();
	for (int row = 0; row<nRows; ++row)
  {
  cout << std::fixed << std::setprecision(3) << inputVector.getElement(row) << endl;
	}
}

template <class T>
int PrintLinEquation(Matrix<T> matrix, Vector<T> vector){
  int numRows = matrix.getNumRows();
  int numCols = matrix.getNumCols();
  if ((numRows != numCols) || (numRows != vector.getNumDim())){
    return 0;
  }
  for (int i = 0; i < numRows; i++){
    for (int j = 0; j < numCols; j++){
      cout << std::fixed << std::setprecision(3) << matrix.GetElement(i, j) << "  ";
    }
    if (i == 1){
      cout << "=  ";
    }else{
      cout << "   ";
    }
    cout << std::fixed << std::setprecision(3) << vector.getElement(i) << endl;
  }
  return 0;
}

int main(){
  cout << "Test sobre resolucion de ecuaciones lineales" << endl;
  cout << "********************************************" << endl;
  cout << endl << "Conversion en matriz escalonada" << endl;
  cout << endl;

  std::vector<double> simpleData = {1.0, 3.0, -1.0, 13.0, 4.0, -1.0, 1.0, 9.0, 2.0, 4.0, 3.0, -6.0};
  Matrix<double> testMatrix(3, 4, &simpleData);

  cout << "Matriz original" << endl;
  testMatrix.printMatrix();
  cout << endl;
  Matrix<double> rowEchelonMatrix = testMatrix.RowEchelon();

  cout << "Matriz escalonada: " << endl;
  rowEchelonMatrix.printMatrix();
  cout << endl;

  std::vector<double> simpleData2 = {1.0, 3.0, -1.0, 4.0, -1.0, 1.0, 2.0, 4.0, 3.0};
  Matrix<double> aMat(3, 3, &simpleData2);
  cout << "Matriz izquierda del sistema de equaciones Ax = b" << endl;
  aMat.printMatrix();
  cout << endl;

  std::vector<double> vectorData= {13.0, 9.0, -6.0};
  Vector<double> bVec(vectorData);
  cout << "Vector derecho del sistema Ax = b" << endl;
  PrintVector(bVec);
  cout << "Sistema en cuestion" << endl;
  PrintLinEquation(aMat, bVec);
  cout << endl;
  int aMatRange = aMat.Rank();
  cout << "El rango de la matriz A es R(A) = " << aMatRange << endl;
  Vector<double> testResult1(3);
  Vector<double> test1 = LinSolve<double>(aMat, bVec, testResult1);
  cout << "La solucion del sistema es: " << endl;
  PrintVector(testResult1);

  std::random_device myRandomDevice;
  std::mt19937 myRandomGenerator(myRandomDevice());
  std::uniform_real_distribution<double> myDistribution(-25.0, 25.0);

  int numUnknowns = 10;

  std::vector<double> coefficientData;
  std::vector<double> unknownData;
  for (int i = 0; i < (numUnknowns * numUnknowns); i++){
    double randomNumber = myDistribution(myRandomGenerator);
    coefficientData.push_back(randomNumber);
  }
  cout << "Matriz random con coeficientes = " << endl;
  Matrix<double> coefficientMatrix(numUnknowns, numUnknowns, &coefficientData);
  coefficientMatrix.printMatrix();

  for (int i = 0; i < (numUnknowns); i++){
    double randomNumber = myDistribution(myRandomGenerator);
    unknownData.push_back(randomNumber);
  }

  cout << "Vector de las soluciones" << endl;
  Vector<double> unknownVector(unknownData);
  PrintVector(unknownVector);
  cout << endl;
  cout << "Multiplicando la matriz por el vector de soluciones obtenemos el vector derecho" << endl;
  Vector<double> unknownResult = coefficientMatrix * unknownVector;
  PrintVector(unknownResult);

  cout << "Resolvemos el sistema..." << endl;
  Vector<double> testResult2(3);
  int test2 = LinSolve<double>(coefficientMatrix, unknownResult, testResult2);
  PrintVector(testResult2);

  cout << "Comprovamos que la resta con nuestro vector de soluciones este cerca de 0" << endl;
  PrintVector(testResult2 - unknownVector);
  cout << endl;

  cout << "Matriz con rango no maximo" << endl;
  std::vector<double> nonInvMatrixData = {3.0, 1.0, 2.0, 5.0, 3.0, 6.0, 2.0, 2.0, 4.0};
  Matrix<double> nonInvMatrix = Matrix(3, 3, &nonInvMatrixData);
  Matrix<double> nonInvEchelon = nonInvMatrix.RowEchelon();
  nonInvMatrix.printMatrix();
  cout << "Calculo del rango Rank = " << nonInvMatrix.Rank() << endl;
  cout << "Matriz con rango no maximo y sin forma escalonada (metodo determinantes)" << endl;
  std::vector<double> nonInvMatrixData2 = {1.0, 0.0, 1.0, 2.0, 0.0, 2.0, 0.0, 1.0, 0.0};
  Matrix<double> nonInvMatrix2 = Matrix(3, 3, &nonInvMatrixData2);
  Matrix<double> nonInvEchelon2 = nonInvMatrix2.RowEchelon();
  nonInvMatrix2.printMatrix();
  cout << "Calculo del rango Rank = " << nonInvMatrix2.Rank() << endl;

  cout << "Resolucion de un sistema de equaciones con infinitas soluciones" << endl;
  cout << "Utilizamos la matriz sin rango maximo con el vector: " << endl;
  std::vector<double> nonInvMatrixData3 = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  Matrix<double> nonInvMatrix3 = Matrix(3, 3, &nonInvMatrixData3);
  std::vector<double> vectorData2 = {1.0, 1.0, 0.0};
  Vector<double> vector2(vectorData2);
  Vector<double> testResult3(3);
  PrintVector(vector2);
  nonInvMatrix2.printMatrix();
  int test3 = LinSolve(nonInvMatrix3, vector2, testResult3);
  cout << "La solucion retornada deberia ser -1: " << test3 << endl;
  vector2.setElement(2, 1.0);
  int test4 = LinSolve(nonInvMatrix3, vector2, testResult3);
  cout << "Si el vector es de solo 1.0 entonces no tiene solucion y el resultado es -2: " << test4 << endl;
  return 0;
}
