
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
#include "../PCA.h"


using namespace std;

int main(){
    cout << "*******************************************" << endl;
    cout << "Test codigo Principal Component Analysis" << endl;
    cout << "*******************************************" << endl;
    cout << endl;

    {
        cout << "Test con 100 observaciones y 3 variables:" << endl;
        string rowData;
        string number;
        stringstream rowDataStream;
        std::vector<double> testData;
        int numRows = 0;
        int numCols = 0;
        ifstream inputFile("Test/PCATestData.csv");

        if (inputFile.is_open()){
            cout << "Fichero abierto correctamente!" << endl;
            while (!inputFile.eof()){
                getline(inputFile, rowData);

                rowDataStream.clear();
                rowDataStream.str(rowData);

                if (numRows < 1)
                    numCols = 0;

                while (rowDataStream.good()){
                    getline(rowDataStream, number, ',');
                    double a = atof(number.c_str());
                    cout << a << endl;
                    testData.push_back(atof(number.c_str()));

                    if (numRows < 1)
                        numCols++;
                }
                numRows++;
            }
        }

        inputFile.close();

        numRows--;
        testData.pop_back();

        cout << "Lectura de fichero completa!" << endl;
        cout << "Leidas " << numRows << " observaciones de " << numCols << " variables." << endl;
        cout << testData.at(0) << " " << testData.at(testData.size() - 1) << endl;
        cout << testData.size() << " elementos en total." << endl;

        Matrix<double> X (numRows, numCols, testData);

        std::vector<double> columnMeans = PCA::ComputeColumnMeans(X);
        Matrix<double> X2 = X;
        PCA::SubractColumnMeans(X2, columnMeans);

        Matrix<double> covX = PCA::ComputeCovariance(X2);
        cout << endl;

        cout << "La matriz de covarianza es: " << endl;
        covX.printMatrix();
        cout << endl;

        Matrix<double> eigenVectors;
        int testResult = PCA::ComputeEigenvector(covX, eigenVectors);
        cout << endl;
        cout << "Los valores propios son: " << endl;
        eigenVectors.printMatrix();

        cout << endl;
        cout << "Testeando la funcion globalmente..." << endl;
        Matrix<double> eigenVectors2;
        int testResult2 = PCA::PCA(X, eigenVectors2);
        cout << "testResult2 = " << testResult2 << endl;
        eigenVectors2.printMatrix();

        cout << endl;
        cout << "Testeando la reduccion dimensional" << endl;
        cout << "La matriz X tiene" << X.getNumRows() << " filas y " << X.getNumCols() << " columnas." << endl;
        cout << endl;
        cout << "Utilizando solo los dos principales componentes:" << endl;
        Matrix<double> V, part2;
        eigenVectors.Separate(V, part2, 2);
        V.printMatrix();
        cout << endl;

        Matrix<double> newX = (V.Transpose() * X.Transpose()).Transpose();
        cout << "El resultado tiene " << newX.getNumRows() << " filas y " << newX.getNumCols() << " columns." << endl;

        ofstream outputFile("Test/PCATestData_Reduced.csv");
        if (outputFile.is_open()){
            for (int i=0; i<newX.getNumRows(); i++){
                outputFile << newX.GetElement(i, 0) << "," << newX.GetElement(i, 1) << endl;
            }
            outputFile.close();
        }
    }
}