#include "Matrix.hpp"

#include <cmath>
#include <cstring>
#include <stdexcept>



// сделать приватными
static bool same_dimension(const Matrix &mat1, const Matrix &mat2) {
    return mat1.rows() == mat2.rows() && mat1.cols() == mat2.cols();
}

static bool is_squared(const Matrix &mat) {
    return mat.rows() == mat.cols();
}

static bool is_compatible(const Matrix &mat1, const Matrix &mat2) {
    return mat1.cols() == mat2.rows();
}
//

static Matrix diagonal(const Matrix &mat, int* sign) {
    if (!is_squared(mat) || !mat.isValid()) {
        return {};
    }

    Matrix diag = mat;

    *sign = -1;
    for (int i = 0; i < diag.rows(); i++) {
        // Поиск главного элемента в столбце
        int maxRow = i;
        for (int k = i + 1; k < diag.rows(); k++) {
            if (fabs(diag.coeffRef(k, i)) > fabs(diag.coeffRef(maxRow,i))) {
                maxRow = k;
            }
        }

        // Меняем строки, если нужно
        if (maxRow != i) {
            for (int k = 0; k < diag.rows(); k++) {
                double temp = diag.coeffRef(i, k);
                diag.coeffRef(i, k) = diag.coeffRef(maxRow, k);
                diag.coeffRef(maxRow, k) = temp;
            }
            *sign = -*sign;
        }

        // Приводим матрицу к верхнетреугольному виду
        for (int k = i + 1; k < diag.rows(); k++) {
            double factor = diag.coeffRef(k, i) / diag.coeffRef(i, i);
            for (int j = i; j < diag.rows(); j++) {
                diag.coeffRef(k, j) -= factor * diag.coeffRef(i, j);
            }
        }
    }
    return diag;
}

// Constructors

Matrix::Matrix()= default;

Matrix::Matrix(size_t cols){
    allocate(1, cols);
}

Matrix::Matrix(size_t rows, size_t cols){
    allocate(rows, cols);
}

Matrix::Matrix(const Matrix &mat){
    if (!mat.isValid()) return;
    allocate(mat.m_rows, mat.m_cols);
    memcpy(m_data, mat.m_data, m_rows*m_cols * sizeof(double));
}

// Destructors

Matrix::~Matrix() {
    deallocate();
}

// Overloaded operators

Matrix Matrix::operator*(const Matrix &mat) const {
    if (!is_compatible(*this, mat) || !isValid() || !mat.isValid()) {
        return {};
    }
    Matrix result(*this);
    result *= mat;

    return result;
}

Matrix Matrix::operator-(const Matrix &mat) const {
    if (!same_dimension(*this, mat) || !isValid() || !mat.isValid()) {
        return {};
    }
    Matrix result(*this);
    result -= mat;
    return result;
}

Matrix Matrix::operator+(const Matrix &mat) const {
    if (!same_dimension(*this, mat) || !isValid() || !mat.isValid()) {
        return {};
    }
    Matrix result(*this);
    result += mat;
    return result;
}

Matrix Matrix::operator*(double value) const {
    Matrix result(*this);

    if (!result.isValid() || !isValid()) return {};

    for (size_t i = 0; i < m_rows * m_cols; ++i) {
        result.m_data[i] *= value;
    }

    return result;
}

Matrix Matrix::operator/(double value) const {
    if (value <= EPS) {
        throw std::runtime_error("division by zero");
    }
    return *this * (1 / value);
}

Matrix & Matrix::operator=(const Matrix &mat) {
    if (mat.isValid()) {
        deallocate();
        allocate(m_rows, m_cols);
        memcpy(m_data, mat.m_data, m_rows*m_cols * sizeof(double));
    } else deallocate();

    return *this;
}

Matrix & Matrix::operator*=(const Matrix &mat) { // наоборот реализовать через *
    if (!is_compatible(*this, mat) || !isValid() || !mat.isValid()) {
        deallocate();
        return *this;
    }

    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < mat.m_cols; ++j) {
            double sum = 0;
            for (size_t k = 0; k < m_cols; ++k) {
                sum += coeffRef(i, k) * mat.coeffRef(k, j);
            }
            coeffRef(i, j) = sum;
        }
    }
    return *this;
}

Matrix & Matrix::operator+=(const Matrix &mat) {
    if (same_dimension(*this, mat) || !isValid() || !mat.isValid()) {
        deallocate();
        return *this;
    }
    for (size_t i = 0; i < m_rows*m_cols; ++i) {
        m_data[i] = m_data[i] + mat.m_data[i];
    }
    return *this;
}

Matrix & Matrix::operator-=(const Matrix &mat) {
    if (same_dimension(*this, mat) || !isValid() || !mat.isValid()) {
        deallocate();
        return *this;
    }
    for (size_t i = 0; i < m_rows*m_cols; ++i) {
        m_data[i] = m_data[i] - mat.m_data[i];
    }
    return *this;
}

Matrix & Matrix::operator*=(const double value) {
    if (!isValid()) {
        return *this;
    }
    *this = *this * value;
    return *this;
}

Matrix & Matrix::operator/=(const double value) {
    if (value <= EPS) {
        deallocate();
        return *this;
    }
    *this *= (1 / value);
    return *this;
}

// Tools

bool Matrix::isValid() const {
    return m_data != nullptr && m_rows > 0 && m_cols > 0;
}

void Matrix::resize(size_t rows, size_t cols) {
    deallocate();
    allocate(rows, cols);
}

// rowIdx * m_cols + colIdx
const double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) const {
    if (rowIdx >= m_rows || colIdx >= m_cols) {
        throw std::out_of_range("out of range");
    }
    if (!isValid()) {
        throw std::out_of_range("invalid matrix data");
    }
    return m_data[rowIdx*m_cols + colIdx];
}

// rowIdx * m_cols + colIdx
double & Matrix::coeffRef(size_t rowIdx, size_t colIdx) {
    if (rowIdx >= m_rows || colIdx >= m_cols) {
        throw std::out_of_range("out of range");
    }
    if (!isValid()) {
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
    for (size_t y = 0; y < m_rows; ++y) {
        for (size_t x = 0; x < m_cols; ++x) {
            coeffRef(y, x) = x == y ? 1.0 : 0.0;
        }
    }
    return *this;
}

Matrix & Matrix::setZero() {
    setConstants(0.0);
    return *this;
}

Matrix & Matrix::setConstants(double value) {
    if (!isValid()) {
        return *this;
    }
    for (size_t i = 0; i < m_rows*m_cols; ++i) {
        m_data[i] = value;
    }
    return *this;
}

Matrix & Matrix::setIdentity(size_t rows, size_t cols) {
    resize(rows, cols);
    setIdentity();
    return *this;
}

Matrix & Matrix::setZero(size_t rows, size_t cols) {
    this->resize(rows, cols);
    setConstants(0);
    return *this;
}

Matrix & Matrix::setConstants(size_t rows, size_t cols, double value) {
    resize(rows, cols);
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
    if (!is_squared(*this) || !isValid() || fabs(det()) < EPS) {
        return {};
    }

    Matrix copy = *this;
    Matrix inverse;
    inverse.resize(m_rows, m_cols);
    for (int i = 0; i < inverse.rows(); i++)
        for (int j = 0; j < inverse.rows(); j++)
            inverse.coeffRef(i, j) = (i == j) ? 1.0 : 0.0;

    for (int i = 0; i < copy.rows(); i++) {
        int maxRow = i;
        for (int k = i + 1; k < copy.rows(); k++)
            if (fabs(copy.coeffRef(k, i)) > fabs(copy.coeffRef(maxRow,i))) {
                maxRow = k;
            }

        // меняем местами строки
        for (int k = 0; k < copy.rows(); k++) {
            double temp = copy.coeffRef(i, k);
            copy.coeffRef(i, k) = copy.coeffRef(maxRow, k);
            copy.coeffRef(maxRow, k) = temp;

            temp = inverse.coeffRef(i, k);
            inverse.coeffRef(i, k) = inverse.coeffRef(maxRow, k);
            inverse.coeffRef(maxRow, k) = temp;
        }

        // Нормализация ведущего элемента
        double diag = copy.coeffRef(i, i);
        for (int k = 0; k < copy.rows(); k++) {
            copy.coeffRef(i, k) /= diag;
            inverse.coeffRef(i, k) /= diag;
        }

        // Обнуление остальных элементов в столбце
        for (int k = 0; k < copy.rows(); k++) {
            if (k == i) continue;
            double factor = copy.coeffRef(k, i);
            for (int j = 0; j < copy.rows(); j++) {
                copy.coeffRef(k, j) -= factor * copy.coeffRef(i, j);
                inverse.coeffRef(k, j) -= factor * inverse.coeffRef(i, j);
            }
        }
    }

    return inverse;
}

double Matrix::det() const{
    if (!is_squared(*this) || !isValid()) {
        return NAN;
    }

    int sign = 1;
    Matrix diag = diagonal(*this, &sign);

    if (!diag.isValid()) return NAN; //

    double det = sign;

    for (int i = 0; i < diag.rows(); i++)
        det *= diag.coeffRef(i,i);

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

void Matrix::allocate(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) {
        return;
    }
    m_cols = cols;
    m_rows = rows;
    m_data = new double[rows * cols];
}

void Matrix::deallocate() {
    m_cols = 0;
    m_rows = 0;
    delete[] m_data;
}
