#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <sstream>
#include <vector>
#include <fstream>

#include "../Matrix.h"

using namespace std;


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
    std::string fileName = "Test/test.csv";
    std::ifstream input(fileName);
    if (input.is_open()){
        cout << "Fichero abierto con exito" << endl;
        std::string currentLine;
        getline(input, currentLine);

        size_t delimiterPosition = currentLine.find(" ");
        std::string nRowString = currentLine.substr(delimiterPosition, currentLine.length());

        int nRows = atoi(nRowString.c_str());
        std::stringstream ss;
        double currentNumber;

        std::vector<double> currentRow;
        std::vector<double> currentMatrixData;
        std::vector<Matrix<double> *> failedTests;

        int totalRows = 0;
        int rowCount = 0;
        int numSuccess = 0;
        int numFails = 0;

        while (!input.eof()){
            getline(input, currentLine);

            ss << currentLine;

            while (ss >> currentNumber){
                currentRow.push_back(currentNumber);
                if (ss.peek() == ',')
                    ss.ignore();
            }

            ss.clear();
            currentMatrixData.insert(currentMatrixData.end(), currentRow.begin(), currentRow.end());
            currentRow.clear();
            rowCount++;
            if (rowCount > 4){
                rowCount = 0;

                Matrix<double> newMatrix(5, 10, currentMatrixData);
                Matrix<double> leftMatrix;
                Matrix<double> rightMhatrix;

                newMatrix.Separate(leftMatrix, rightMhatrix, 5);

                leftMatrix.Inverse();
                if (leftMatrix.Compare(rightMhatrix, 1e-9)){
                    numSuccess++;
                }else{
                    cout << "Invert test failed" << endl;
                    cout << numSuccess << endl;
                    numFails++;
                    failedTests.push_back(new Matrix<double>(newMatrix));
                }
                currentMatrixData.clear();
            }
            totalRows++;
        }
        cout << "Num succes: " << numSuccess << endl;
        cout << "Num fails: " << numFails << endl;
    }
    input.close();
    return 0;
}