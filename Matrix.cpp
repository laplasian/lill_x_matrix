#include "Matrix.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

// Constructors:

Matrix::Matrix(){
}

Matrix::Matrix(size_t cols): m_rows(1), m_cols(cols){
    m_data = allocate(1, cols);
}

Matrix::Matrix(size_t rows, size_t cols): m_rows(rows), m_cols(cols){
    m_data = allocate(m_rows, m_cols);
}

Matrix::Matrix(int rows, int cols, const std::vector<double> &values): m_rows(rows), m_cols(cols){
    m_data = allocate(m_rows, m_cols);
    if (!isValid() || values.size() != m_rows * m_cols) return;
    memcpy(m_data, values.data(), sizeof(double) * m_rows * m_cols);
}

Matrix::Matrix(const Matrix &mat): m_rows(mat.m_rows), m_cols(mat.m_cols) {
    m_data = allocate(m_rows, m_cols);
    if (!isValid() || !mat.isValid() || m_rows != mat.m_rows || m_cols != mat.m_cols) return;
    memcpy(m_data, mat.m_data, m_rows*m_cols * sizeof(double));
}

// Destructors

Matrix::~Matrix() {
    m_data = deallocate(m_data);
}

// Overloaded operators

Matrix Matrix::operator*(const Matrix &mat) const {
    if (m_cols != mat.m_rows || !isValid() || !mat.isValid()) {
        return {};
    }
    Matrix result(m_rows, mat.m_cols);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < mat.m_cols; ++j) {
            double sum = 0;
            for (size_t k = 0; k < m_cols; ++k) {
                sum += coeffRef(i, k) * coeffRef(k, j);
            }
            result.coeffRef(i, j) = sum;
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
    result.m_data = new double[m_rows * m_cols];
    for (size_t i = 0; i < m_rows * m_cols; ++i) {
        result.m_data[i] *= value;
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
        delete m_data;
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
    *this = *this * mat;
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

// rowIdx * m_cols + colIdx
const double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) const {
    if (rowIdx >= m_rows || colIdx >= m_cols) {
        throw std::out_of_range("out of range");
    } else if (m_isValid == false || m_data == nullptr) {
        throw std::out_of_range("invalid matrix data");
    }
    return m_data[rowIdx*m_cols + colIdx];
}


// rowIdx * m_cols + colIdx
double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) {
    if (rowIdx >= m_rows || colIdx >= m_cols) {
        throw std::out_of_range("out of range");
    } else if (m_isValid == false || m_data == nullptr) {
        throw std::out_of_range("invalid matrix data");
    }
    return m_data[rowIdx*m_cols + colIdx];
}

double & Matrix::operator()(size_t rowIdx, size_t colIdx) {
    return coeffRef(rowIdx, colIdx);
}

const double & Matrix::operator()(size_t rowIdx, size_t colIdx) const {
    return coeffRef(rowIdx, colIdx);
}

bool Matrix::operator==(const Matrix &mat) const {
    if (m_cols != mat.m_cols || m_rows != mat.m_rows || m_isValid == false || mat.m_isValid == false) {
        return false;
    }
    for (size_t i = 0; i < m_cols * m_rows; ++i) {
        if (m_data[i] != mat.m_data[i]) return false;
    }
    return true;
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

Matrix Matrix::transpose() const{
    Matrix _new(m_cols, m_rows);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            _new.coeffRef(j, i) = coeffRef(i, j);
        }
    }
    return _new;
}

Matrix Matrix::inverse() const{
    Matrix _new(m_rows, m_cols);
    if (m_rows != m_cols || m_isValid == false || det() == 0 || m_data == nullptr) {
        return {};
    }

    return _new;
}

double Matrix::det() const{
    if (m_rows != m_cols || m_isValid == false || m_data == nullptr) {
        return NAN;
    }

    Matrix copy = *this;

    double det = 1;
    int sign = 1;  // Считаем количество перестановок строк

    for (int i = 0; i < m_rows; i++) {
        // Поиск главного элемента в столбце
        int maxRow = i;
        for (int k = i + 1; k < m_rows; k++) {
            if (fabs(copy.coeffRef(k, i)) > fabs(copy.coeffRef(maxRow,i))) {
                maxRow = k;
            }
        }

        // Если главный элемент - ноль, определитель равен 0
        if (fabs(copy.coeffRef(maxRow, i)) < 1e-9) {
            return 0.0;
        }

        // Меняем строки, если нужно
        if (maxRow != i) {
            for (int k = 0; k < m_rows; k++) {
                double temp = coeffRef(i, k);
                copy.coeffRef(i, k) = copy.coeffRef(maxRow, k);
                copy.coeffRef(maxRow, k) = temp;
            }
            sign = -sign;  // Меняем знак определителя при перестановке строк
        }

        // Приводим матрицу к верхнетреугольному виду
        for (int k = i + 1; k < m_rows; k++) {
            double factor = copy.coeffRef(k, i) / copy.coeffRef(i, i);
            for (int j = i; j < m_rows; j++) {
                copy.coeffRef(k, j) -= factor * copy.coeffRef(i, j);
            }
        }

        // Умножаем диагональный элемент в определитель
        det *= copy.coeffRef(i, i);
    }

    return det * sign;

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

double * Matrix::allocate(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) {
        return nullptr;
    }
    return new double[rows*cols];
}

double * Matrix::deallocate(double *data) {
    if (data == nullptr) {
        return nullptr;
    }
    delete[] data;
    return nullptr;
}