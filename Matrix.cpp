#include "Matrix.h"

#include <string.h>

// Constructors:

Matrix::Matrix(): m_rows(0), m_cols(0), m_data(nullptr), m_isValid(false) {
}

Matrix::Matrix(size_t cols): m_rows(1), m_cols(cols), m_data(nullptr), m_isValid(true) {
    if (cols == 0) {
        m_isValid = false;
    }
    m_data = new double[cols];
}

Matrix::Matrix(size_t rows, size_t cols): m_rows(rows), m_cols(cols), m_data(nullptr), m_isValid(true) {
    if (cols * rows == 0) {
        m_isValid = false;
    }
    m_data = new double[m_rows*m_cols];
}

Matrix::Matrix(const Matrix &mat): m_rows(mat.m_rows), m_cols(mat.m_cols), m_data(nullptr) {
    if (mat.m_data == nullptr || mat.m_rows * mat.m_cols == 0 || mat.m_isValid == false) {
        m_isValid = false;
        return;
    }
    m_data = new double[m_rows*m_cols];
    memcpy(m_data, mat.m_data, m_rows*m_cols * sizeof(double));
}

// Destructors

Matrix::~Matrix() {
    delete[] m_data;
}

// Overwrite operators

Matrix Matrix::operator*(const Matrix &mat) const {
    if (m_cols != mat.m_rows || m_isValid == false || mat.m_isValid == false || m_data == nullptr || mat.m_data == nullptr) {
        return {};
    }
    Matrix result(m_rows, mat.m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < mat.m_cols; ++j) {
            double sum = 0;
            for (size_t k = 0; k < m_cols; ++k) {
                sum += m_data[i*m_cols + k] * mat.m_data[k*m_cols + j];
            }
            result.m_data[i*m_cols + j] = sum;
        }
    }

    return result;
}

Matrix Matrix::operator-(const Matrix &mat) const {
    return *this + mat * (-1);
}

Matrix Matrix::operator+(const Matrix &mat) const {
    if (m_cols != mat.m_cols || m_rows != mat.m_rows || m_isValid == false || mat.m_isValid == false || m_data == nullptr || mat.m_data == nullptr) {
        return {};
    }
    Matrix result(m_rows, mat.m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < mat.m_cols; ++j) {
            result.m_data[(i*m_cols + j)] = m_data[i*m_cols + j] + mat.m_data[i*m_cols + j];
        }
    }
    return result;
}

Matrix Matrix::operator*(double value) const {
    if (m_isValid == false || m_data == nullptr) {
        return {};
    }
    Matrix result(m_rows, m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            result.m_data[i*m_cols + j] = m_data[i*m_cols + j] * value;
        }
    }
    return result;
}

Matrix Matrix::operator/(double value) const {
    return *this * (1 / value);
}

Matrix & Matrix::operator=(const Matrix &mat) {
    m_rows = mat.m_rows;
    m_cols = mat.m_cols;
    m_isValid = mat.m_isValid;
    if (mat.m_data != nullptr) {
        delete[] m_data;
        m_data = new double[m_rows*m_cols];
        memcpy(m_data, mat.m_data, m_rows*m_cols * sizeof(double));
    } else m_data = nullptr;

    return *this;
}

Matrix & Matrix::operator*=(const Matrix &mat) {
}

Matrix & Matrix::operator+=(const Matrix &mat) {
}

Matrix & Matrix::operator-=(const Matrix &mat) {
}

Matrix & Matrix::operator*=(double value) {
}

Matrix & Matrix::operator/=(double value) {
}

// Tools

bool Matrix::isValid() {
}

void Matrix::resize(size_t rows, size_t cols) {
}

const double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) const {
}

double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) {
}

const double * Matrix::data() const {
}

double * Matrix::data() {
}

size_t Matrix::rows() {
}

size_t Matrix::cols() {
}

Matrix & Matrix::setIdentity() {
}

Matrix & Matrix::setZero() {
}

Matrix & Matrix::setConstants(double value) {
}

Matrix & Matrix::setIdentity(size_t rows, size_t cols) {
}

Matrix & Matrix::setZero(size_t rows, size_t cols) {
}

Matrix & Matrix::setConstants(size_t rows, size_t cols, double value) {
}

Matrix Matrix::transpose() {
}

Matrix Matrix::inverse() {
}

double Matrix::det() {
}

Matrix Matrix::identity(size_t rows, size_t cols) {
}

Matrix Matrix::zeros(size_t rows, size_t cols) {
}

Matrix Matrix::constants(size_t rows, size_t cols, double value) {
}

Matrix operator*(double value, const Matrix &mat) {
    return value *= mat;
}
