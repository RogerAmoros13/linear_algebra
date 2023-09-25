#ifndef LINSOLVE_H
#define LINSOLVE_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "Matrix.h"
#include "Vector.h"

using namespace std;

constexpr int LINSOLVE_NOUNIQUESOLUTION = -1;
constexpr int LINSOLVE_NOSOLUTIONS = -2;

template <class T>
int LinSolve(const Matrix<T> &aMatrix, const Vector<T> &bVector, Vector<T> &resultVec){
    Matrix<T> inputMatrix = aMatrix;

    int originalRank = inputMatrix.Rank();

    int numDims = bVector.getNumDim();
    std::vector<T> vecData;
    for (int i = 0; i < numDims; i++){
        vecData.push_back(bVector.getElement(i));
    }

    Matrix<T> bMatrix(numDims, 1, &vecData);

    inputMatrix.join(bMatrix);

    int augmentedRank = inputMatrix.Rank();
    if ((originalRank == augmentedRank) && (originalRank < inputMatrix.getNumRows())){
        return LINSOLVE_NOUNIQUESOLUTION;
    }else if (originalRank < augmentedRank){
        return LINSOLVE_NOSOLUTIONS;
    }else{
        Matrix<T> rowEchelonMatrix = inputMatrix.RowEchelon();
        Vector<T> output(vecData);

        int numRows = rowEchelonMatrix.getNumRows();
        int numCols = rowEchelonMatrix.getNumCols();
        int startRow = numRows - 1;
        for (int i = startRow; i >= 0; i--){
            T currentResult = rowEchelonMatrix.GetElement(i, numCols - 1);

            T cumSum = static_cast<T>(0.0);
            for (int j = i + 1; j < numRows; j++){
                cumSum += (rowEchelonMatrix.GetElement(i, j) * output.getElement(j));
            }

            T result = ( currentResult - cumSum) / rowEchelonMatrix.GetElement(i, i);

            output.setElement(i, result);
        }
        resultVec = output;
    }
    return 1;
}

#endif
