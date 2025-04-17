#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <vector>

#define EPS 1e-6

class Matrix final {
public:
    // Constructors
    Matrix();

    explicit Matrix(size_t cols);
    Matrix(size_t rows, size_t cols);
    ~Matrix();

    Matrix(const Matrix& mat);

    Matrix operator*(const Matrix& mat) const;
    Matrix operator-(const Matrix& mat) const;
    Matrix operator+(const Matrix& mat) const;

    Matrix operator*(double value) const;
    Matrix operator/(double value) const;

    Matrix& operator=(const Matrix& mat);
    Matrix& operator*=(const Matrix& mat);
    Matrix& operator+=(const Matrix& mat);
    Matrix& operator-=(const Matrix& mat);

    Matrix& operator*=(double value);
    Matrix& operator/=(double value);

    bool isValid() const;

    void resize(size_t rows, size_t cols);

    const double& coeffRef(size_t rowIdx, size_t colIdx) const;
    double& coeffRef(size_t rowIdx, size_t colIdx);

    const double* data() const;
    double* data();

    size_t rows() const;
    size_t cols() const;

    Matrix& setIdentity();
    Matrix& setZero();
    Matrix& setConstants(double value);

    Matrix& setIdentity(size_t rows, size_t cols);
    Matrix& setZero(size_t rows, size_t cols);
    Matrix& setConstants(size_t rows, size_t cols, double value);

    Matrix transpose() const;
    Matrix inverse() const;
    double det() const;

    static Matrix identity(size_t rows, size_t cols);
    static Matrix zeros(size_t rows, size_t cols);
    static Matrix constants(size_t rows, size_t cols, double value);

    friend Matrix operator*(double value, const Matrix& mat);
private:
    size_t m_rows{};
    size_t m_cols{};
    double* m_data{};

    bool same_dimension(const Matrix &mat) const;
    bool is_squared() const;
    bool is_compatible(const Matrix &mat) const;

    static Matrix diagonal(const Matrix &mat, int* sign);

    void allocate(size_t rows, size_t cols);
    void deallocate();
};

#endif //MATRIX_H
