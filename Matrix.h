#ifndef MATRIX_H
#define MATRIX_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "Vector.h"

template <class T>
class Matrix{
    public:
        Matrix();
        Matrix(int nRows, int nCols);
        Matrix(const Matrix<T> &inputMatrix);
        Matrix(int nRows, int nCols, const T *inputData);
        Matrix(int nRows, int nCols, const std::vector<T> &inputData);

        // Destructor
        ~Matrix();

        // Configuration Method
        bool resize(int numRows, int numCols);
        void setToIdendity();

        // Funciones para obtener y escribir un elemento de la matriz.
        T GetElement(int row, int col) const;
        bool setElement(int row, int col, T elementValue);

        // Funciones para obtener el número de columnas/filas.
        int getNumRows() const;
        int getNumCols() const;

        // Compute matrix inverse

        bool Inverse();
        Matrix<T> RowEchelon();
        T Determinant();
        int Rank();
        Matrix<T> Transpose() const;


        // Sobreescribir los operadores para trabajar con mayor facilidad.
        bool operator == (const Matrix<T>& rElement);
        bool Compare (const Matrix<T>& matrix1, double epsilon);

        Matrix<T> operator = (const Matrix<T>& lElement);

        template <class U> friend Matrix<U> operator + (const Matrix<U>& lElement, const Matrix<U>& rElement);
        template <class U> friend Matrix<U> operator + (const Matrix<U>& lElement, const U& rElement);
        template <class U> friend Matrix<U> operator + (const U& lElement, const Matrix<U>& rElement);

        template <class U> friend Matrix<U> operator - (const Matrix<U>& lElement, const Matrix<U>& rElement);
        template <class U> friend Matrix<U> operator - (const Matrix<U>& lElement, const U& rElement);
        template <class U> friend Matrix<U> operator - (const U& lElement, const Matrix<U>& rElement);

        template <class U> friend Matrix<U> operator * (const Matrix<U>& lElement, const Matrix<U>& rElement);
        template <class U> friend Matrix<U> operator * (const Matrix<U>& lElement, const U& rElement);
        template <class U> friend Matrix<U> operator * (const U& lElement, const Matrix<U>& rElement);

        template <class U> friend Vector<U> operator * (const Matrix<U>& lElement, const Vector<U>& rElement);

        bool Separate(Matrix<T> &matrix1, Matrix<T> &matrix2, int nCol);

        static int Rank(const Matrix<T> &matrix);

    // public:
    // private:
        int Sub2Ind(int row, int col) const;
        bool isSquare();
        bool isSymmetric();
        bool isRowEchelon();
        bool isNonZero();
        bool closeEnough(T value1, T value2);
        void swapRow(int row1, int row2);
        void multRow(int row, T factor);
        void multAdd(int row1, int row2, T factor);
        bool join(const Matrix<T>& matrix2);
        int findRowWithMaxElement(int colNumber, int startingRow);
        void printMatrix();
        void printMatrix(int rows);
        Matrix<T> FindSubMatrix(int rowNum, int colNum);

    private:
        T *m_matrixData;
        int m_nRows, m_nCols, m_nElements;
};

// DEFINICIÓN DE LAS FUNCIONES DECLARADAS

// Constructor

// Basic constructor
template <class T>
Matrix<T>::Matrix(){
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++){
        m_matrixData[i] = 0.0;
    }
}

// Construct zero matriz with nRows and nCols
template <class T>
Matrix<T>::Matrix(int nRows, int nCols){
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = nRows * nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++){
        m_matrixData[i] = 0.0;
    }
}

// Construct a matrix from another matrix
template <class T>
Matrix<T>::Matrix(const Matrix<T> &inputMatrix)
{
	m_nRows = inputMatrix.m_nRows;
	m_nCols = inputMatrix.m_nCols;
	m_nElements = inputMatrix.m_nElements;

	m_matrixData = new T[m_nElements];
	for (int i=0; i<m_nElements; i++)
		m_matrixData[i] = inputMatrix.m_matrixData[i];
}

// Construct a matrix giving nRows, nCols and the data.
template <class T>
Matrix<T>::Matrix(int nRows, int nCols, const T *inputData){
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i){
        m_matrixData[i] = inputData[i];
    }
}

// Construct matrix from std::vector
template <class T>
Matrix<T>::Matrix(int nRows, int nCols, const std::vector<T> &inputData){
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = nRows * nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++){
        m_matrixData[i] = inputData.at(i);
    }
}
// Destructor

template <class T>
Matrix<T>::~Matrix()
{
	// Destructor.
	if (m_matrixData)
		delete[] m_matrixData;
	m_matrixData = nullptr;
}

// Funciones

template <class T>
bool Matrix<T>::resize(int numRows, int numCols){
    m_nRows = numRows;
    m_nCols = numCols;
    m_nElements = numCols * numRows;
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    if (m_matrixData != nullptr){
        for (int i = 0; i < m_nElements; i++){
            m_matrixData[i] = 0.0;
        }
        return true;
    }
    return false;
}

template <class T>
T Matrix<T>::GetElement(int row, int col) const{
    int pos = Sub2Ind(row, col);
    if (pos >= 0){
        return m_matrixData[pos];
    }
    return 0.0;
}

template <class T>
bool Matrix<T>::setElement(int row, int col, T elementValue){
    int index = Sub2Ind(row, col);
    if ((index >= 0) && (index < m_nElements)){
        m_matrixData[index] = elementValue;
        return true;
    }
    return false;
}

template <class T>
int Matrix<T>::getNumCols() const{
    return m_nCols;
}


template <class T>
int Matrix<T>::getNumRows() const{
    return m_nRows;
}

// **************************************************************************************
// Funciones de los operadores

// Operador comparativo
template <class T>
bool Matrix<T>::operator== (const Matrix<T>& rElement){
    if ((this->m_nRows != rElement.m_nRows) && (this->m_nCols != rElement.m_nCols)){
        return false;
    }
    for (int i = 0; i < this->m_nElements; i++){
        // if (this->m_matrixData[i] != rElement.m_matrixData[i]){
        if (!closeEnough(this->m_matrixData[i], rElement.m_matrixData[i])){
            return false;
        }
    }
    return true;
}

template <class T>
Matrix<T> Matrix<T>::operator = (const Matrix<T> &rhs){
    if (this != &rhs){
        m_nRows = rhs.m_nRows;
        m_nCols = rhs.m_nCols;
        m_nElements = rhs.m_nElements;
        if (m_matrixData){
            delete[] m_matrixData;
        }
        m_matrixData = new T[m_nElements];
        for (int i = 0; i < m_nElements; i++){
            m_matrixData[i] = rhs.m_matrixData[i];
        }
    }
    return *this;
}

template <class T>
Vector<T> operator * (const Matrix<T>& lElement, const Vector<T>& rElement){
    if (lElement.m_nCols != rElement.getNumDim()){
        throw std::invalid_argument("Las dimensiones de la matriz y el vector no coinciden");
    }

    Vector<T> result(lElement.m_nRows);
    for (int i = 0; i < lElement.m_nRows; i++){
        T cumSum = static_cast<T>(0.0);
        for (int j = 0; j < lElement.m_nCols; j++){
            cumSum += lElement.GetElement(i, j) * rElement.getElement(j);
        }
        result.setElement(i, cumSum);
    }
    return result;
}

// Suma matriz + matriz
template <class T>
Matrix<T> operator + (const Matrix<T>& lElement, const Matrix<T>& rElement){
    int numRows = lElement.m_nRows;
    int numCols = lElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement.m_matrixData[i] + rElement.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Suma constante + Matriz
template <class T>
Matrix<T> operator + (const T& lElement, const Matrix<T>& rElement){
    int numRows = rElement.m_nRows;
    int numCols = rElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement + rElement.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Suma Matriz + constante
template <class T>
Matrix<T> operator + (const Matrix<T>& lElement, const T& rElement){
    int numRows = lElement.m_nRows;
    int numCols = lElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement.m_matrixData[i] + rElement;
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Resta matriz - matriz
template <class T>
Matrix<T> operator - (const Matrix<T>& lElement, const Matrix<T>& rElement){
    int numRows = lElement.m_nRows;
    int numCols = lElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement.m_matrixData[i] - rElement.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Resta constante - Matriz
template <class T>
Matrix<T> operator - (const T& lElement, const Matrix<T>& rElement){
    int numRows = rElement.m_nRows;
    int numCols = rElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement - rElement.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Resta Matriz - constante
template <class T>
Matrix<T> operator - (const Matrix<T>& lElement, const T& rElement){
    int numRows = lElement.m_nRows;
    int numCols = lElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement.m_matrixData[i] - rElement;
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Multiplicación Matriz * Matriz
template <class T>
Matrix<T> operator * (const Matrix<T>& lElement, const Matrix<T>& rElement){
    int n = lElement.m_nRows;
    int m = lElement.m_nCols;
    int m2 = rElement.m_nRows;
    int p = rElement.m_nCols;
    if (m == m2){
        T *tempResult = new T[lElement.m_nRows * rElement.m_nCols];
        for (int i = 0; i < n; i++){
            for (int j = 0; j < p; j++){
                T elementResult = 0.0;
                for (int k = 0; k < m; k++){
                    int lIndex = i * m + k;
                    int rIndex = k * p + j;
                    elementResult += (lElement.m_matrixData[lIndex] * rElement.m_matrixData[rIndex]);
                }
                tempResult[i * p + j] = elementResult;
            }
        }
        Matrix<T> result(n, p, tempResult);
        delete[] tempResult;
        return result;
    }else{
        Matrix<T> result(1, 1);
        return result;
    }
}

// Multiplicación const * matriz
template <class T>
Matrix<T> operator * (const T& lElement, const Matrix<T>& rElement){
    int numRows = rElement.m_nRows;
    int numCols = rElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement * rElement.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Multiplicación Matriz * constante
template <class T>
Matrix<T> operator * (const Matrix<T>& lElement, const T& rElement){
    int numRows = lElement.m_nRows;
    int numCols = lElement.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++){
        tempResult[i] = lElement.m_matrixData[i] * rElement;
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}


// *********************************************************************************
// Funciones privadas - publicas

template <class T>
int Matrix<T>::Sub2Ind(int row, int col) const{
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0)){
        return row * m_nCols + col;
    }
    return -1;
}

template <class T>
void Matrix<T>::setToIdendity(){
    if (!isSquare()){
        throw std::invalid_argument("No se puede formar una identidad si la matriz no es cuadada");
    }
    for (int i = 0; i < m_nRows; i++){
        for (int j = 0; j < m_nCols; j++){
            if (i == j){
                m_matrixData[Sub2Ind(i, j)] = 1.0;
            }
            else{
                m_matrixData[Sub2Ind(i, j)] = 0.0;
            }
        }
    }
}

template <class T>
bool Matrix<T>::Compare (const Matrix<T>& matrix1, double epsilon){
    int numRows1 = matrix1.m_nRows;
    int numCols1 = matrix1.m_nCols;
    if ((numRows1 != m_nRows) || (numCols1 != m_nCols)){
        return false;
    }
    double cumSum = 0.0;
    for (int i = 0; i < m_nElements; i++){
        T a = m_matrixData[i];
        T b = matrix1.m_matrixData[i];
        cumSum += ((a - b) * (a - b));
    }
    double finalValue = sqrt(cumSum / ((numRows1 * numCols1) - 1));
    if (finalValue < epsilon){
        return true;
    }
    return false;

}
template <class T>
bool Matrix<T>::isSquare(){
    if (m_nRows == m_nCols){
        return true;
    }else{
        return false;
    }
}

template <class T>
bool Matrix<T>::isSymmetric(){
    if (!this->isSquare()){
        return false;
    }
    T currentRowElement = static_cast<T>(0.0);
    T currentColElement = static_cast<T>(0.0);
    bool returnFlag = true;
    int diagIndex = 0;

    while ((diagIndex < m_nCols) && returnFlag){
        int rowIndex = diagIndex + 1;
        while ((rowIndex < m_nRows) && returnFlag){
            currentRowElement = this->GetElement(rowIndex, diagIndex);
            currentColElement = this->GetElement(diagIndex, rowIndex);

            if (!closeEnough(currentRowElement, currentColElement)){
                returnFlag = false;
            }
            rowIndex++;
        }
        diagIndex++;
    }
    return returnFlag;
}

template <class T>
bool Matrix<T>::closeEnough(T value1, T value2){
    return fabs(value1 - value2) < 1e-9;
}

template <class T>
void Matrix<T>::swapRow(int row1, int row2){
    T *tempRow = new T[m_nCols];
    for (int i = 0; i < m_nCols; i++){
        tempRow[i] = m_matrixData[Sub2Ind(row1, i)];
    }
    for (int i = 0; i < m_nCols; i++){
        m_matrixData[Sub2Ind(row1, i)] = m_matrixData[Sub2Ind(row2, i)];
    }
    for (int i = 0; i < m_nCols; i++){
        m_matrixData[Sub2Ind(row2, i)] = tempRow[i];
    }
    delete[] tempRow;
}

template <class T>
void Matrix<T>::multRow(int row, T factor){
    for (int i = 0; i < m_nCols; i++){
        m_matrixData[Sub2Ind(row, i)] = m_matrixData[Sub2Ind(row, i)] * factor;
    }
}

template <class T>
void Matrix<T>::multAdd(int row1, int row2, T factor){
    for (int i = 0; i < m_nCols; i++){
        m_matrixData[Sub2Ind(row1, i)] += m_matrixData[Sub2Ind(row2, i)] * factor;
    }
}

template <class T>
int Matrix<T>::findRowWithMaxElement(int colNumber, int startingRow){
    T tempValue = m_matrixData[Sub2Ind(startingRow, colNumber)];
    int rowIndex = startingRow;
    for (int i = startingRow + 1; i < m_nRows; i++){
        if (fabs(m_matrixData[Sub2Ind(i, colNumber)]) > fabs(tempValue)){
            tempValue = m_matrixData[Sub2Ind(i, colNumber)];
            rowIndex = i;
        }
    }
    return rowIndex;
}

template <class T>
bool Matrix<T>::join(const Matrix<T>& matrix2){
    int numRows1 = m_nRows;
    int numRows2 = matrix2.m_nRows;
    int numCols1 = m_nCols;
    int numCols2 = matrix2.m_nCols;

    if (numRows1 != numRows2){
        throw std::invalid_argument("No se pueden juntar matrices con numero distinto de filas");
    }

    T* newMatrixData = new T[numRows1 * (numCols1 + numCols2)];
    int linearIndex, resultLinearIndex;
    for (int i = 0; i < numRows1; i++){
        for (int j = 0; j < (numCols1 + numCols2); j++){
            resultLinearIndex = i * (numCols1 + numCols2) + j;

            if (j < numCols1){
                linearIndex = i * numCols1 + j;
                newMatrixData[resultLinearIndex] = m_matrixData[linearIndex];
            }else{
                linearIndex = i * numCols2 + j - numCols1;
                newMatrixData[resultLinearIndex] = matrix2.m_matrixData[linearIndex];
            }
        }
    }
    m_nCols = numCols1 + numCols2;
    m_nElements = numRows1 * (numCols1 + numCols2);
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++){
        m_matrixData[i] = newMatrixData[i];
    }
    delete[] newMatrixData;
    return true;
}

template <class T>
bool Matrix<T>::Inverse(){
    if (!isSquare()){
        throw std::invalid_argument("No es una matriz cuadrada");
    }
    Matrix<T> identityMatrix(m_nRows, m_nCols);
    identityMatrix.setToIdendity();

    int originalNumCols = m_nCols;
    join(identityMatrix);

    int cRow, cCol;
    int maxCount = 100;
    int count = 0;
    bool completeFlag = false;
    while ((!completeFlag) && (count < maxCount)){
        for (int diagIndex = 0; diagIndex < m_nRows; diagIndex++){
            cRow = diagIndex;
            cCol = diagIndex;

            int maxIndex = findRowWithMaxElement(cCol, cRow);
            if (maxIndex != cRow){
                swapRow(cRow, maxIndex);
            }

            if (m_matrixData[Sub2Ind(cRow, cCol)] != 1.0){
                T multFactor = 1.0 / m_matrixData[Sub2Ind(cRow, cCol)];
                multRow(cRow, multFactor);
            }

            for (int rowIndex = cRow + 1; rowIndex < m_nRows; rowIndex++){
                if (!closeEnough(m_matrixData[Sub2Ind(rowIndex, cCol)], 0.0)){
                    int rowOneIndex = cCol;

                    T currentElementValue = m_matrixData[Sub2Ind(rowIndex, cCol)];
                    T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, cCol)];

                    if (!closeEnough(rowOneValue, 0.0)){
                        T correctionFactor = -(currentElementValue / rowOneValue);
                        multAdd(rowIndex, rowOneIndex, correctionFactor);
                    }
                }
            }
            for (int colIndex = cCol + 1; colIndex < originalNumCols; colIndex++){
                if (!closeEnough(m_matrixData[Sub2Ind(cRow, colIndex)], 0.0)){
                    int rowOneIndex = colIndex;

                    T currentElementValue = m_matrixData[Sub2Ind(cRow, colIndex)];
                    T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, colIndex)];

                    if (!closeEnough(rowOneValue, 0.0)){
                        T correctionFactor = -(currentElementValue / rowOneValue);
                        multAdd(cRow, rowOneIndex, correctionFactor);
                    }
                }
            }
        }
        Matrix<T> leftHalf;
        Matrix<T> rightHalf;
        this->Separate(leftHalf, rightHalf, originalNumCols);
        // leftHalf.printMatrix();
        if (leftHalf == identityMatrix){
            completeFlag = true;
            m_nCols = originalNumCols;
            m_nElements = m_nRows * m_nCols;
            delete[] m_matrixData;
            m_matrixData = new T[m_nElements];
            for (int i = 0; i < m_nElements; i++){
                m_matrixData[i] = rightHalf.m_matrixData[i];
            }

        }
        count++;
    }
    return completeFlag;
}

template <class T>
bool Matrix<T>::Separate(Matrix<T> &matrix1, Matrix<T> &matrix2, int nCol){
    int numRows = m_nRows;
    int numCols1 = nCol;
    int numCols2 = m_nCols - nCol;

    matrix1.resize(numRows, numCols1);
    matrix2.resize(numRows, numCols2);

    for (int i = 0; i < m_nRows; i++){
        for (int j = 0; j < m_nCols; j++){
            if (j < nCol){
                matrix1.setElement(i, j, this->GetElement(i, j));
            }
            else{
                matrix2.setElement(i, j - nCol, this->GetElement(i, j));
            }
        }
    }
    return true;
}

template <class T>
T Matrix<T>::Determinant(){
    if (!isSquare()){
        throw std::invalid_argument("No es cuadrada!");
    }
    if (m_nElements == 4){
        return m_matrixData[0] * m_matrixData[3] - m_matrixData[1] * m_matrixData[2];
    }else{
        int sign = 1;
        double cumulativeSum = 0;
        for (int j = 0; j < m_nCols; j++){
            Matrix<T> subMatrix = this->FindSubMatrix(0, j);
            cumulativeSum += this->GetElement(0, j) * subMatrix.Determinant() * sign;
            sign = -sign;
        }
        return cumulativeSum;
    }
}

template <class T>
Matrix<T> Matrix<T>::FindSubMatrix(int rowNum, int colNum){
    Matrix<T> subMatrix(m_nRows - 1, m_nCols - 1);
    int count = 0;
    for (int i = 0; i < m_nRows; i++){
        for (int j = 0; j < m_nCols; j++){
            if ((i != rowNum) && (j != colNum)){
                subMatrix.m_matrixData[count] = this->GetElement(i, j);
                count++;
            }
        }
    }
    return subMatrix;
}

template <class T>
bool Matrix<T>::isRowEchelon(){
    T cumSum = static_cast<T>(0.0);
    for (int i = 0; i < m_nRows; i++){
        for (int j = 0; j < i; j++){
            cumSum += m_matrixData[Sub2Ind(i, j)];
        }
    }
    return closeEnough(cumSum, 0.0);
}

template <class T>
bool Matrix<T>::isNonZero(){
    T cumSum = static_cast<T>(0.0);
    bool nonZero = false;
    for (int i = 0; i < m_nElements; i++){
        if (!closeEnough(m_matrixData[i], 0.0)){
            nonZero = true;
        }
    }
    return nonZero;
}

template <class T>
Matrix<T> Matrix<T>::RowEchelon(){
    if (m_nCols < m_nRows){
        throw std::invalid_argument("La matriz tiene que tener igual o más columnas que filas.");
    }

    T *tempMatrixData;
    tempMatrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++){
        tempMatrixData[i] = m_matrixData[i];
    }
    int cRow, cCol;
    int maxCount = 100;
    int count = 0;
    bool completeFlag = false;
    while ((count < maxCount) && (!completeFlag)){
        for (int diagIndex = 0; diagIndex < m_nRows; diagIndex++){
            cRow = diagIndex;
            cCol = diagIndex;

            for (int rowIndex = cRow+1; rowIndex < m_nRows; rowIndex++){
                if (!closeEnough(m_matrixData[Sub2Ind(rowIndex, cCol)], 0.0)){
                    int rowOneIndex = cCol;
                    T currentElementValue = m_matrixData[Sub2Ind(rowIndex, cCol)];
                    T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, cCol)];
                    if (!closeEnough(rowOneValue, 0.0)){
                        T correctionFactor = -(currentElementValue / rowOneValue);
                        multAdd(rowIndex, rowOneIndex, correctionFactor);
                    }
                }
            }
        }
        completeFlag = this->isRowEchelon();
        count++;
    }
    Matrix<T> outputMatrix(m_nRows, m_nCols, m_matrixData);
    for (int i = 0; i < m_nElements; i++){
        m_matrixData[i] = tempMatrixData[i];
    }
    return outputMatrix;
}

template <class T>
int Matrix<T>::Rank(){
    Matrix<T> matrixCopy = this->RowEchelon();
    int numNonZeroRows = 0;
    if (!matrixCopy.isRowEchelon()){
        std::vector<Matrix<T>> subMatrixVector;
        subMatrixVector.push_back(*this);
        bool completeFlag = false;
        int subMatrixCount = 0;
        while ((subMatrixCount < subMatrixVector.size()) && (!completeFlag)){
            Matrix<T> currentMatrix = subMatrixVector[subMatrixCount];
            subMatrixCount++;

            if (currentMatrix.isNonZero()){
                T currentMatrixDet = currentMatrix.Determinant();
                if (!closeEnough(currentMatrixDet, 0.0)){
                    completeFlag = true;
                    numNonZeroRows = currentMatrix.getNumRows();
                }else{
                    if ((currentMatrix.getNumCols() > 2) && (currentMatrix.getNumRows() > 2)){
                        for (int i = 0; i < currentMatrix.getNumCols(); i++){
                            for (int j = 0; j < currentMatrix.getNumRows(); j++){
                                subMatrixVector.push_back(currentMatrix.FindSubMatrix(i, j));
                            }
                        }
                    }
                }
            }
        }
    }else{
        int nRows = matrixCopy.getNumRows();
        int nCols = matrixCopy.getNumCols();

        for (int i = 0; i < nRows;i++){
            int colSum = 0.0;
            for (int j = 0; j < nCols; j++){
                if (!closeEnough(matrixCopy.GetElement(i, j), 0.0)){
                    colSum++;
                }
            }
            if (colSum > 0){
                numNonZeroRows++;
            }
        }
    }
    return numNonZeroRows;
}

template <class T>
Matrix<T> Matrix<T>::Transpose() const{
    Matrix<T> resultMatrix(m_nCols, m_nRows);
    for (int i = 0; i < m_nRows; i++){
        for (int j = 0; j < m_nCols; j++){
            resultMatrix.setElement(j, i, this->m_matrixData[Sub2Ind(i, j)]);
        }
    }
    return resultMatrix;
}

template <class T>
void Matrix<T>::printMatrix(){
    int rows = getNumRows();
    int cols = getNumCols();
    for (int i = 0; i < rows; i++){
        std::cout << endl;
        for (int j = 0; j < cols; j++){
            std::cout << std::fixed << std::setprecision(3) << GetElement(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

#endif
