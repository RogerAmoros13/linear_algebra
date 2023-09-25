#ifndef LSQ_H
#define LSQ_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "Matrix.h"
#include "Vector.h"

constexpr int LSQ_NOINVERSE = -1;

template <class T>
int LSQ(const Matrix<T> &Xin, const Vector<T> &yin, Vector<T> &result){
    Matrix<T> X = Xin;
    Vector<T> y = yin;

    Matrix<T> XT = X.Transpose();

    Matrix<T> XXT = XT * X;

    if (!XXT.Inverse()){
        return LSQ_NOINVERSE;
    }

    Matrix<T> XXTXT = XXT * XT;

    result = XXTXT * y;
    return 1;
}

#endif
