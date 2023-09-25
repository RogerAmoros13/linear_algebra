#ifndef EIG_H
#define EIG_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "Matrix.h"
#include "Vector.h"
#include "QR.h"

constexpr int EIG_MATRIXNOTSQUARE = -1;
constexpr int EIG_MAXITERATIONSEXCEEDED = -2;
constexpr int EIG_MATRIXNOTSYMMETRIC = -3;

template <typename T>
int EigQR(const Matrix<T> &inputMatrix, std::vector<T> &eigenValues){
    Matrix<T> A = inputMatrix;

    if (!A.isSquare()){
        return EIG_MATRIXNOTSQUARE;
    }

    if (!A.isSymmetric()){
        return EIG_MATRIXNOTSYMMETRIC;
    }

    int numRows = A.getNumRows();

    Matrix<T> identityMatrix(numRows, numRows);
    identityMatrix.setToIdendity();

    Matrix<T> Q (numRows, numRows);
    Matrix<T> R (numRows, numRows);

    int maxIterations = 10e3;
    int iterationCount = 0;
    bool continueFlag = true;

    while ((iterationCount < maxIterations) && continueFlag){
        int returnValue = QR<T>(A, Q, R);

        A = R * Q;

        if (A.isRowEchelon()){
            continueFlag = false;
        }
        iterationCount++;
    }

    for (int i=0; i<numRows; i++){
        eigenValues.push_back(A.GetElement(i, i));
    }

    if (iterationCount == maxIterations){
        return EIG_MAXITERATIONSEXCEEDED;
    }else{
        return 0;
    }
}

template <typename T>
int InvPIt(const Matrix<T> &inputMatrix, const T &eigenValue, Vector<T> &eigenVector){
    Matrix<T> A = inputMatrix;

    if (!A.isSquare()){
        return EIG_MATRIXNOTSQUARE;
    }

    std::random_device myRandomDevice;
    std::mt19937 myRandomGenerator(myRandomDevice());
    std::uniform_int_distribution<int> myDistribution(1.0, 10.0);

    int numRows = A.getNumRows();

    Matrix<T> identityMatrix(numRows, numRows);
    identityMatrix.setToIdendity();

    Vector<T> v(numRows);
    for (int i=0; i<numRows; i++){
        v.setElement(i, static_cast<T>(myDistribution(myRandomGenerator)));
    }

    int maxIterations = 100;
    int iterationCount = 0;
    T deltaThreshold = static_cast<T>(1e-9);
    T delta = static_cast<T>(1e6);
    Vector<T> prevVector(numRows);
    Matrix<T> tempMatrix(numRows, numRows);

    while ((iterationCount < maxIterations) && (delta > deltaThreshold)){
        prevVector = v;
        tempMatrix = A - (eigenValue*identityMatrix);
        tempMatrix.printMatrix();
        tempMatrix.Inverse();
        tempMatrix.printMatrix();
        v = tempMatrix * v;
        v.Normalize();

        delta = (v - prevVector).norm();

        iterationCount++;
    }

    eigenVector = v;

    if (iterationCount == maxIterations){
        return EIG_MAXITERATIONSEXCEEDED;
    }else{
        return 0;
    }
}

template <class T>
int EIG_PIt(const Matrix<T> X, T &eigenValue, Vector<T> &eigenVector){
    Matrix<T> inputMatrix = X;

    if (!inputMatrix.isSquare()){
        return EIG_MATRIXNOTSQUARE;
    }

    std::random_device myRandomDevice;
    std::mt19937 myRandomGenerator(myRandomDevice());
    std::uniform_int_distribution<int> myDistribution(1.0, 10.0);

    int numRows = inputMatrix.getNumRows();

    Matrix<T> identityMatrix(numRows, numRows);
    identityMatrix.setToIdendity();

    Vector<T> v(numRows);
    for (int i = 0; i < numRows; i++){
        v.setElement(i, static_cast<T>(myDistribution(myRandomGenerator)));
    }

    Vector<T> v1(numRows);
    int numIterations = 1000;
    for (int i = 0; i < numIterations; i++){
        v1 = inputMatrix * v;
        v1.Normalize();
        v = v1;
    }

    eigenVector = v1;
    T cumSum = static_cast<T>(0.0);
    for (int i = 1; i < numRows; i++){
        cumSum += inputMatrix.GetElement(0, i) * v1.getElement(i);
    }
    eigenValue = (cumSum / v1.getElement(0)) + inputMatrix.GetElement(0, 0);

    return 0;
}

#endif
