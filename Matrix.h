#ifndef MATRIX_H
#define MATRIX_H

class Matrix final {
public:
    // Constructors
    Matrix();
    Matrix(size_t cols);
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

    bool isValid();

    void resize(size_t rows, size_t cols);

    const double& coeffRef(size_t rowIdx, size_t colIdx) const;
    double& coeffRef(size_t rowIdx, size_t colIdx);

    const double* data() const;
    double* data();

    size_t rows();
    size_t cols();

    Matrix& setIdentity();
    Matrix& setZero();
    Matrix& setConstants(double value);

    Matrix& setIdentity(size_t rows, size_t cols);
    Matrix& setZero(size_t rows, size_t cols);
    Matrix& setConstants(size_t rows, size_t cols, double value);

    Matrix transpose();
    Matrix inverse();
    double det();

    static Matrix identity(size_t rows, size_t cols);
    static Matrix zeros(size_t rows, size_t cols);
    static Matrix constants(size_t rows, size_t cols, double value);

    friend Matrix operator*(double value, const Matrix& mat);

private:
    bool m_isValid{};
    size_t m_rows{};
    size_t m_cols{};
    double* m_data{};
};

#endif //MATRIX_H
