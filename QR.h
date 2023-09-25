#ifndef QR_H
#define QR_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "Matrix.h"
#include "Vector.h"

using namespace std;

constexpr int QR_MATRIXNOTSQUARE = -1;

template <typename T>
int QR(const Matrix<T> &A, Matrix<T> &Q, Matrix<T> &R)
{
    Matrix<T> inputMatrix = A;

    if (!inputMatrix.isSquare())
    {
        return QR_MATRIXNOTSQUARE;
    }

    int numCols = inputMatrix.getNumCols();
    std::vector<Matrix<T>> PList;

    for (int j = 0; j < (numCols - 1); j++)
    {
        Vector<T> a1(numCols - j);
        Vector<T> b1(numCols - j);

        for (int i = j; i < numCols; i++)
        {
            a1.setElement(i - j, inputMatrix.GetElement(i, j));
            b1.setElement(i - j, static_cast<T>(0.0));
        }
        b1.setElement(0, static_cast<T>(1.0));

        T a1norm = a1.norm();

        int sgn = -1;
        if (a1.getElement(0) < static_cast<T>(0.0))
        {
            sgn = 1;
        }

        Vector<T> u = a1 - (sgn * a1norm * b1);
        Vector<T> n = u.Normalized();

        Matrix<T> nMat(numCols - j, 1);
        for (int i = 0; i < (numCols - j); i++)
        {
            nMat.setElement(i, 0, n.getElement(i));
        }

        Matrix<T> nMatT = nMat.Transpose();
        Matrix<T> I(numCols - j, numCols - j);
        I.setToIdendity();

        Matrix<T> Ptemp = I - (static_cast<T>(2.0) * nMat) * nMatT;

        Matrix<T> P(numCols, numCols);
        P.setToIdendity();
        for (int row = j; row < numCols; row++)
        {
            for (int col = j; col < numCols; col++)
            {
                P.setElement(row, col, Ptemp.GetElement(row - j, col - j));
            }
        }

        PList.push_back(P);
        inputMatrix = P * inputMatrix;

        int size = PList.size();
    }

    Matrix<T> Qmat = PList.at(0);
    for (int i = 1; i < (numCols - 1); i++)
    {
        Qmat = Qmat * PList.at(i).Transpose();
    }

    Q = Qmat;

    int numElements = PList.size();
    Matrix<T> Rmat = PList.at(numElements - 1);
    for (int i = (numElements - 2); i >= 0; i--)
    {
        Rmat = Rmat * PList.at(i);
    }
    Rmat = Rmat * A;

    R = Rmat;
    return 1;
}

#endif
