#ifndef VECTOR_H
#define VECTOR_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>


template <class T>
class Vector{
    public:
        Vector();
        Vector(std::vector<T> inputData);
        Vector(int numDims);
        ~Vector();


        int getNumDim() const;
        T getElement(int index) const;
        void setElement(int index, T value);

        T norm();

        Vector<T> Normalized();
        void Normalize();

        // Operacions
        Vector<T> operator + (const Vector<T> &rhs) const;
        Vector<T> operator - (const Vector<T> &rhs) const;
        Vector<T> operator * (const T &rhs) const;

        template <class U> friend Vector<U> operator* (const U &lhs, const Vector<U> &rhs);

        static T dot(const Vector<T> &a, const Vector<T> &b);
        static Vector<T> cross(const Vector<T> &a, const Vector<T> &b);

        void PrintVector();
    private:
        std::vector<T> m_vectorData;
        int m_nDims;
};


template <class T>
void Vector<T>::PrintVector()
{
	int dims = getNumDim();
	for (int i = 0; i < dims; ++i){
    std::cout << std::fixed << std::setprecision(3) << m_vectorData[i] << std::endl;
	}
    std::cout << std::endl;
}

template <class T>
Vector<T>::Vector(){
    m_nDims = 0;
    m_vectorData = std::vector<T>();
}

template <class T>
Vector<T>::Vector(std::vector<T> inputData){
    m_nDims = inputData.size();
    m_vectorData = inputData;
}

template <class T>
Vector<T>::Vector(int numDims){
    m_nDims = numDims;
    m_vectorData = std::vector<T>(numDims, static_cast<T>(0.0));
}

template <class T>
Vector<T>::~Vector(){}

template <class T>
int Vector<T>::getNumDim() const{
    return m_nDims;
}

template <class T>
T Vector<T>::getElement(int index) const{
    return m_vectorData[index];
}

template <class T>
void Vector<T>::setElement(int index, T value){
    m_vectorData[index] = value;
}

template <class T>
Vector<T> Vector<T>::operator+ (const Vector<T> &rhs) const{
    if (m_nDims != rhs.m_nDims){
        throw std::invalid_argument("Dimensiones distintas");
    }

    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; i++){
        resultData.push_back(m_vectorData.at(i) + rhs.m_vectorData.at(i));
    }
    Vector<T> result(resultData);
    return result;
}

template <class T>
Vector<T> Vector<T>::operator- (const Vector<T> &rhs) const{
    if (m_nDims != rhs.m_nDims){
        throw std::invalid_argument("Dimensiones distintas");
    }

    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; i++){
        resultData.push_back(m_vectorData.at(i) - rhs.m_vectorData.at(i));
    }
    Vector<T> result(resultData);
    return result;
}

template <class T>
Vector<T> Vector<T>::operator* (const T &rhs) const{
    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; i++){
        resultData.push_back(m_vectorData.at(i) * rhs);
    }
    Vector<T> result(resultData);
    return result;
}

template <class T>
Vector<T> operator* (const T &lhs, const Vector<T> &rhs){
    std::vector<T> resultData;
    for (int i = 0; i < rhs.m_nDims; i++){
        resultData.push_back(lhs * rhs.m_vectorData.at(i));
    }
    Vector<T> result(resultData);
    return result;
}

template <class T>
T Vector<T>::dot(const Vector<T> &a, const Vector<T> &b){
    if (a.m_nDims != b.m_nDims){
        throw std::invalid_argument("Dimensiones distintas!");
    }
    T cumulativeSum = 0;
    for (int i = 0; i < a.m_nDims; i++){
        cumulativeSum += a.m_vectorData[i] * b.m_vectorData[i];
    }
    return cumulativeSum;
}

template <class T>
Vector<T> Vector<T>::cross(const Vector<T> &a, const Vector<T> &b){
    if (a.m_nDims != b.m_nDims){
        throw std::invalid_argument("Dimensiones distintas!");
    }

    if (a.m_nDims != 3){
        throw std::invalid_argument("Dimensiones distintas!");
    }

    std::vector<T> resultData;
    resultData.push_back(a.m_vectorData.at(1) * b.m_vectorData.at(2) - a.m_vectorData.at(2) * b.m_vectorData.at(1));
    resultData.push_back(a.m_vectorData.at(0) * b.m_vectorData.at(2) - a.m_vectorData.at(2) * b.m_vectorData.at(0));
    resultData.push_back(a.m_vectorData.at(0) * b.m_vectorData.at(1) - a.m_vectorData.at(1) * b.m_vectorData.at(0));

    Vector<T> result(resultData);
    return result;
}

template <class T>
T Vector<T>::norm(){
    T cumSum = static_cast<T>(0.0);
    for (int i = 0; i < m_nDims; i++){
        cumSum += (m_vectorData.at(i) * m_vectorData.at(i));
    }
    return sqrt(cumSum);
}

template <class T>
Vector<T> Vector<T>::Normalized(){
    T vecNorm = this->norm();
    Vector<T> result(m_vectorData);
    return result * (static_cast<T>(1.0) / vecNorm);
}

template <class T>
void Vector<T>::Normalize(){
    T vecNorm = this->norm();
    for (int i = 0; i < m_nDims; i++){
        m_vectorData.at(i) = m_vectorData.at(i) / vecNorm;
    }
}

#endif
