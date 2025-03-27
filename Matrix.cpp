#include "Matrix.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

// Constructors:

Matrix::Matrix(): m_rows(0), m_cols(0), m_data(nullptr), m_isValid(false) {
}

Matrix::Matrix(size_t cols): m_rows(1), m_cols(cols), m_data(nullptr), m_isValid(true) {
    if (cols == 0) {
        m_isValid = false;
        return;
    }
    m_data = new double[cols];
}

Matrix::Matrix(size_t rows, size_t cols): m_rows(rows), m_cols(cols), m_data(nullptr), m_isValid(true) {
    if (cols * rows == 0) {
        m_isValid = false;
        return;
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
    if (m_cols != mat.m_rows || m_isValid == false || mat.m_isValid == false || m_data == nullptr || mat.m_data == nullptr) {
        m_isValid = false;
        return *this;
    }
    const Matrix result = *this * mat;
    *this = result;
    return *this;
}

Matrix & Matrix::operator+=(const Matrix &mat) {
    if (m_cols != mat.m_cols || m_rows != mat.m_rows || m_isValid == false || mat.m_isValid == false || m_data == nullptr || mat.m_data == nullptr) {
        m_isValid = false;
        return *this;
    }
    const Matrix result = *this + mat;
    *this = result;
    return *this;
}

Matrix & Matrix::operator-=(const Matrix &mat) {
    *this = *this - mat;
    return *this;
}

Matrix & Matrix::operator*=(double value) {
    if (m_isValid == false || m_data == nullptr) {
        return *this;
    }
    const Matrix result = *this * value;
    *this = result;
    return *this;
}

Matrix & Matrix::operator/=(double value) {
    if (value == 0) {
        m_isValid = false;
        return *this;
    }
    *this *= (1 / value);
    return *this;
}

// Tools

bool Matrix::isValid() const {
    return m_isValid;
}

void Matrix::resize(size_t rows, size_t cols) {
    *this = Matrix(rows, cols);
}

const double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) const {
    if (rowIdx >= m_rows || colIdx >= m_cols) {
        throw std::out_of_range("out of range");
    } else if (m_isValid == false || m_data == nullptr) {
        throw std::out_of_range("invalid matrix data");
    }
    return m_data[rowIdx*m_cols + colIdx];
}

double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) {
    if (rowIdx >= m_rows || colIdx >= m_cols) {
        throw std::out_of_range("out of range");
    } else if (m_isValid == false || m_data == nullptr) {
        throw std::out_of_range("invalid matrix data");
    }
    return m_data[rowIdx*m_cols + colIdx];
}

const double * Matrix::data() const {
    return m_data;
}

double * Matrix::data() {
    return m_data;
}

size_t Matrix::rows() const {
    return m_rows;
}

size_t Matrix::cols() const {
    return m_cols;
}

Matrix & Matrix::setIdentity() {
    setConstants(1.0);
    return *this;
}

Matrix & Matrix::setZero() {
    setConstants(0.0);
    return *this;
}

Matrix & Matrix::setConstants(double value) {
    if (m_data == nullptr) {
        return *this;
    }
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            m_data[i*m_cols + j] = value;
        }
    }
    m_isValid = true;
    return *this;
}

Matrix & Matrix::setIdentity(size_t rows, size_t cols) {
    this->resize(rows, cols);
    setConstants(1.0);
    return *this;
}

Matrix & Matrix::setZero(size_t rows, size_t cols) {
    this->resize(rows, cols);
    setConstants(0);
    return *this;
}

Matrix & Matrix::setConstants(size_t rows, size_t cols, double value) {
    this->resize(rows, cols);
    setConstants(value);
    return *this;
}

Matrix Matrix::transpose() {
    Matrix old = *this;
    this->resize(m_cols, m_rows);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            m_data[i*m_cols + j] = old.m_data[j*m_rows + i];
        }
    }
    return *this;
}

Matrix Matrix::inverse() {
    return {};
}

double Matrix::det() {
    if (rows() != cols()) return NAN; // Детерминант определён только для квадратных матриц

    size_t n = rows();
    Matrix temp(*this); // Создаём копию матрицы, чтобы не изменять оригинал
    double det = 1.0;

    for (size_t i = 0; i < n; i++) {
        size_t pivot = i;

        // Ищем строку с максимальным элементом в i-м столбце
        for (size_t j = i + 1; j < n; j++) {
            if (fabs(temp.coeffRef(j, i)) > fabs(temp.coeffRef(pivot, i))) {
                pivot = j;
            }
        }

        // Если ведущий элемент ноль → детерминант ноль
        if (fabs(temp.coeffRef(pivot, i)) < 1e-9) return 0.0;

        // Меняем строки местами (если нужно) и меняем знак детерминанта
        if (i != pivot) {
            for (size_t k = 0; k < n; k++) {
                std::swap(temp.coeffRef(i, k), temp.coeffRef(pivot, k));
            }
            det *= -1;
        }

        // Умножаем детерминант на ведущий элемент
        det *= temp.coeffRef(i, i);

        // Приведение столбца ниже главной диагонали к нулям
        for (size_t j = i + 1; j < n; j++) {
            double coef = temp.coeffRef(j, i) / temp.coeffRef(i, i);
            for (size_t k = i; k < n; k++) {
                temp.coeffRef(j, k) -= coef * temp.coeffRef(i, k);
            }
        }
    }
    return det;
}

Matrix Matrix::identity(size_t rows, size_t cols) {
    Matrix result(rows, cols);
    result.setIdentity();
    return result;
}

Matrix Matrix::zeros(size_t rows, size_t cols) {
    Matrix result(rows, cols);
    result.setZero();
    return result;
}

Matrix Matrix::constants(size_t rows, size_t cols, double value) {
    Matrix result(rows, cols);
    result.setConstants(value);
    return result;
}

Matrix operator*(double value, const Matrix &mat) {
    Matrix result(mat);
    result *= value;
    return result;
}
