#ifndef PCA_H
#define PCA_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>

#include "Matrix.h"
#include "Vector.h"
#include "EIG.h"

constexpr int PCA_MATRIXNOTSQUARE = -1;
constexpr int PCA_MATRIXNOTSYMMETRIC = -2;

namespace PCA{
    template <typename T>
    std::vector<T> ComputeColumnMeans(const Matrix<T> &inputData){
        int numRows = inputData.getNumRows();
        int numCols = inputData.getNumCols();

        std::vector<T> output;

        for (int j=0; j<numCols; j++){
            T cumSum = static_cast<T>(0.0);
            for (int i=0; i<numRows; i++){
                cumSum += inputData.GetElement(i, j);
            }
        output.push_back(cumSum / static_cast<T>(numRows));
        }
        return output;
    }

    template <typename T>
    void SubractColumnMeans(Matrix<T> &inputData, std::vector<T> &columnMeans){
        int numRows = inputData.getNumRows();
        int numCols = inputData.getNumCols();

        for (int j=0; j<numCols; j++){
            for (int i=0; i<numRows; i++){
                inputData.setElement(i, j, inputData.GetElement(i, j) - columnMeans.at(j));
            }
        }
    }

    template <typename T>
    Matrix<T> ComputeCovariance(const Matrix<T> &X){
        int numRows = X.getNumRows();
        Matrix<T> covX = (
            static_cast<T>(1.0) / static_cast<T>(numRows - 1)
        ) * (X.Transpose() * X);
        return covX;
    }

    template <typename T>
    int ComputeEigenvector(const Matrix<T> &covarianceMatrix, Matrix<T> &eigenVectors){
        Matrix<T> X = covarianceMatrix;

        if (!X.isSquare())
            return PCA_MATRIXNOTSQUARE;

        if (!X.isSymmetric())
            return PCA_MATRIXNOTSYMMETRIC;

        std::vector<T> eigenValues;
        int returnStatus = EigQR(X, eigenValues);

        std::sort(eigenValues.begin(), eigenValues.end());
        std::reverse(eigenValues.begin(), eigenValues.end());

        Vector<T> eV(X.getNumCols());
        Matrix<T> eVM(X.getNumRows(), X.getNumCols());
        for (int j=0; j<eigenValues.size(); j++){
            T eig = eigenValues.at(j);
            int returnStatus2 = InvPIt<T>(X, eig, eV);
            for (int i=0; i<eV.getNumDim(); i++){
                eVM.setElement(i, j, eV.getElement(i));
            }
        }

        eigenVectors = eVM;

        return returnStatus;
    }

    template <typename T>
    int PCA(const Matrix<T> &inputData, Matrix<T> &outputComponents){
        Matrix<T> X = inputData;
        std::vector<T> columnMeans = ComputeColumnMeans(X);

        SubractColumnMeans<T>(X, columnMeans);

        Matrix<T> covX = ComputeCovariance(X);

        Matrix<T> eigenVectors;

        int returnStatus = ComputeEigenvector(covX, eigenVectors);

        outputComponents = eigenVectors;

        return returnStatus;
    }
}



#endif
