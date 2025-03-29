#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <vector>

#define EPS 1e-6

class Matrix final {
public:
    // Constructors
    Matrix();
    Matrix(size_t cols);
    Matrix(size_t rows, size_t cols);
    Matrix(int rows, int cols, const std::vector<double>& values); // Дополнение
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

    double& operator()(size_t rowIdx, size_t colIdx); // Дополнение
    const double& operator()(size_t rowIdx, size_t colIdx) const; // Дополнение
   // const double& operator[](size_t rowIdx, size_t colIdx)  const; // Дополнение
    bool operator==(const Matrix& mat) const; // Дополнение


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

    static double *allocate(size_t rows, size_t cols);
    static double * deallocate(double *data);
};

#endif //MATRIX_H
