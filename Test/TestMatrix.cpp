#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>

#include "../Matrix.h"
#include "../Vector.h"

using namespace std;

template <class T>
void PrintVector(Vector<T> vector)
{
	int dims = vector.getNumDim();
	for (int i = 0; i < dims; ++i){
    cout << std::fixed << std::setprecision(3) << vector.getElement(i) << endl;
	}
    cout << endl;
}

template <class T>
void PrintMatrix(Matrix<T> matrix){
    int rows = matrix.getNumRows();
    int cols = matrix.getNumCols();
    for (int i = 0; i < rows; i++){
        cout << "                         ";
        for (int j = 0; j < cols; j++){
            cout << setprecision(3) << matrix.GetElement(i, j) << " ";
        }
        cout << endl;
    }
}


int main(){
    cout << "Codigo para probar las funcionalidades basicas de la clase Matrix" << endl;
    double simpleData[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    Matrix<double> matrix(3, 4, simpleData);
    cout << "***********************  Matriz 3x4  ****************************" << endl;
    PrintMatrix(matrix);

    cout << endl << "*********************  Test atacar a elemento matriz ********************";
    cout << "Entrada (0, 0) = " << matrix.GetElement(0, 0) << endl;
    cout << "Entrada (1, 0) = " << matrix.GetElement(1, 0) << endl;
    cout << "Entrada (0, 1) = " << matrix.GetElement(0, 1) << endl;
    cout << "Entrada (2, 2) = " << matrix.GetElement(2, 2) << endl;
    cout << "Entrada (5, 5) = " << matrix.GetElement(5, 5) << endl;

    cout << endl << "***********************  Multiplicacion  ********************" << endl;
    double simpleData2[12] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    Matrix<double> matrix2(4, 3, simpleData2);
    cout << "***********************  Matriz2 4x3  ****************************" << endl;
    PrintMatrix(matrix2);
    cout << "Matriz * Matriz2 = " << endl;
    Matrix m_m2_mult = matrix * matrix2;
    PrintMatrix(m_m2_mult);

    cout << "Test Multiplication Vector * Matrix" << endl;

    double vectorData[3] = {1.5, 10.0, 3.0};
    Matrix<double> vector(3, 1, vectorData);
    cout << "***********************  Matriz 3x1  ****************************" << endl;
    PrintMatrix(vector);

    double identityData[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    Matrix<double> identity3(3, 3, identityData);
    cout << "***********************  Matriz identidad (3)  ****************************" << endl;
    PrintMatrix(identity3);
    cout << endl << "Resultado Vector * Identidad(3)" << endl;
    PrintMatrix(identity3 * vector);
    cout << endl << "Resultado de Matriz * Escalar(2)" << endl;
    PrintMatrix(matrix * 2.0);
    cout << endl << "Resultado de Escalar(3) * Matriz" << endl;
    PrintMatrix(2.0 * matrix);
    cout << endl << "Resultado Matriz + Matriz" << endl;
    PrintMatrix(matrix + matrix);
    cout << endl << "Resultado Matriz + 4" << endl;
    PrintMatrix(matrix + 4.0);
    cout << endl << "Resultado 4 + Matriz" << endl;
    PrintMatrix(4.0 + matrix);
    cout << endl << "Resultado Matriz - Matriz" << endl;
    PrintMatrix(matrix - matrix);
    cout << endl << "Resultado Matriz - 2.5" << endl;
    PrintMatrix(matrix - 2.5);
    cout << endl << "Resultado 1.2 - Matriz" << endl;
    PrintMatrix(1.2 - matrix);
    cout << endl << "Comprobacion I(3) es cuadrada" << endl;
    cout << identity3.isSquare() << endl;
    double squareData[16] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
    Matrix<double> squareMatrix(4, 4, squareData);
    Matrix<double> squareMatrix2(4, 4, squareData);
    cout << endl << "Matriz 4x4 para ver permutaciones" << endl;
    squareMatrix.printMatrix();
    cout << endl << "Permutacion de fila 1-3 de m(4x4)" << endl;
    squareMatrix.swapRow(0, 2);
    squareMatrix.printMatrix();
    cout << endl << "Multiplicar fila 3 por 2 de m(4x4)" << endl;
    squareMatrix.multRow(2, 3);
    squareMatrix.printMatrix();
    cout << endl << "Sumar 2x fila 2 a la fila 4" << endl;
    squareMatrix.multAdd(3, 1, 2);
    squareMatrix.printMatrix();
    cout << endl << "Buscar el elemento mas alto de la fila 3" << endl;
    cout << squareMatrix.findRowWithMaxElement(2, 0) << endl;
    cout << endl << "Convertir a matriz1 identidad" << endl;
    squareMatrix2.setToIdendity();
    squareMatrix2.printMatrix();
    // squareMatrix.join(squareMatrix2);
    // squareMatrix.printMatrix();
    cout << endl << "Swap a la matriu identitat" << endl;
    identity3.swapRow(0, 2);
    identity3.printMatrix();
    cout << endl << "Matriz 4x4" << endl;
    squareMatrix.printMatrix();
    cout << endl << "Matriz inversa de matriz 4x4" << endl;
    squareMatrix.Inverse();
    squareMatrix.printMatrix();
    double invertTestData2[25] =
			{2.0, -3.0, 4.0, 5.0, 6.0,
			 1.0, 2.0, 3.0, 4.0, 5.0,
			 9.0, 5.0, 3.0, 2.0, 6.0,
			 2.0, 4.0, 6.0, 5.0, 1.0,
			 1.0, 7.0, 5.0, 2.0, 3.0};
    Matrix<double> invertTest = Matrix(5, 5, invertTestData2);
    Matrix<double> invertTest2 = invertTest;
    cout << endl << " lskjdfanlskfjanljk nl" << endl;
    (invertTest2*invertTest).printMatrix();
    bool flag = invertTest.Inverse();
    cout << "************ " << flag << " *************" << endl;
    invertTest.printMatrix();
    cout << endl;
    (invertTest2*invertTest).printMatrix();
    cout << endl << "Prova definitiva " << endl;
    double data3inverse[9] = {
        2.0, 1.0, 1.0,
        1.0, 0.0, 1.0,
        0.0, 3.0, 1.0
    };
    Matrix<double> matrix3inverse = Matrix(3, 3, data3inverse);
    matrix3inverse.printMatrix();
    cout << endl << "La inversa es: " << endl;
    matrix3inverse.Inverse();
    matrix3inverse.printMatrix();

    cout << "Mi rollito" << endl;

    std::vector<double> my_data = {0.264, 0.005, 0.393, 0.005, 0.249, 0.352, 0.393, 0.352, 1.062};
    Matrix<double> my_matrix(3, 3, my_data);
    my_matrix.Inverse();
    my_matrix.printMatrix();
    return 0;
}