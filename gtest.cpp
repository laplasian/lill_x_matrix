#include <gtest/gtest.h>
#include "Matrix.hpp"

void set(Matrix& mat, std::vector<double> v) {
    if (v.size() != mat.cols()*mat.rows()) return;
    memcpy(mat.data(), v.data(), sizeof(double)*mat.cols()*mat.rows());
};

bool operator==(const Matrix &a, const Matrix &b) {
    if (!a.isValid() || !b.isValid() || a.cols()!= b.cols() || a.rows() != b.rows()) return false;
    for (int i = 0; i < a.cols()*a.rows(); i++) {
        if (a.data()[i] != b.data()[i]) return false;
    }
    return true;
}


// TEST CONSTRUCTORS

TEST(MatrixTest, Constructors) {
    Matrix m1(2, 2);
    EXPECT_EQ(m1.rows(), 2);
    EXPECT_EQ(m1.cols(), 2);
    set(m1, {1,2,3,4});
    // m1.coeffRef(0,0) = 1.0;
    // m1.coeffRef(0, 1) = 2.0;
    // m1.coeffRef(1, 0) = 3.0;
    // m1.coeffRef(1, 1) = 4.0;
    // EXPECT_DOUBLE_EQ(m1.coeffRef(0, 0), 1.0);
    // EXPECT_DOUBLE_EQ(m1.coeffRef(0, 1), 2.0);
    // EXPECT_DOUBLE_EQ(m1.coeffRef(1, 0), 3.0);
    // EXPECT_DOUBLE_EQ(m1.coeffRef(1, 1), 4.0);

    Matrix m2(2);
    EXPECT_EQ(m2.rows(), 1);
    EXPECT_EQ(m2.cols(), 2);

    Matrix m3(m1);
    EXPECT_EQ(m3.rows(), m1.rows());
    EXPECT_EQ(m3.cols(), m1.cols());
    EXPECT_EQ(m1, m3);
}

// TEST BaseOperations

TEST(MatrixTest, BaseOperations) {
    Matrix m1(2, 2);
    Matrix m2(2, 2);
    Matrix m3(2, 2);
    Matrix m4(2, 2);
    set(m1, {1,0,1,0});
    set(m2, {1,0,1,0});
    set(m3, {0,0,0,0});
    set(m4, {2,0,2,0});

    EXPECT_EQ(m1 * m2, m1);
    EXPECT_EQ(m1 - m2, m3);
    EXPECT_EQ(m1 + m2, m4);
    EXPECT_EQ(m1 * 2, m4);
    EXPECT_EQ(m4 / 2, m1);
}

// TEST Операции с присваиванием

TEST(MatrixTest, ReturnOperations) {
    Matrix m1(2, 2);
    Matrix m2(2, 2);
    Matrix m3(2, 2);
    Matrix m4(2, 2);
    set(m1, {1,0,1,0});
    set(m2, {1,0,1,0});
    set(m3, {0,0,0,0});
    set(m4, {2,0,2,0});

    Matrix result(2,2);

    result = m1;
    EXPECT_EQ(result, m1);
    result += m2;
    EXPECT_EQ(result, m4);
    result -= m3;
    EXPECT_EQ(result, m4);
    result *= m1;
    EXPECT_EQ(result, m4);
    result /= 2;
    EXPECT_EQ(result, m1);
    result *= 2;
    EXPECT_EQ(result, m4);
}

// TEST TOOLS

TEST(MatrixTest, TOOLS) {
    Matrix m1(2, 2);
    EXPECT_EQ(m1.rows(), 2);
    EXPECT_EQ(m1.cols(), 2);
    set(m1, {1,2,3,4});
    m1.coeffRef(0,0) = 1.0;
    m1.coeffRef(0, 1) = 2.0;
    m1.coeffRef(1, 0) = 3.0;
    m1.coeffRef(1, 1) = 4.0;
    EXPECT_DOUBLE_EQ(m1.coeffRef(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m1.coeffRef(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(m1.coeffRef(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(m1.coeffRef(1, 1), 4.0);

    std::vector<double> scr = {1, 2, 3, 4};
    Matrix m2(2,2);
    set(m2, scr);
    EXPECT_EQ(scr.data(), m2.data());

    m1.resize(3,3);
    EXPECT_EQ(m1.rows(), 3);
    EXPECT_EQ(m1.cols(), 3);
}

// TEST SET FUNCTIONS

TEST(MatrixTest, SetFunctions) {
    Matrix m1(2, 2);
    Matrix m2(2, 2);

    m1.setIdentity();
    set(m2, {1,0,1,0});
    EXPECT_EQ(m1, m2);

    m1.setZero();
    set(m2, {0,0,0,0});
    EXPECT_EQ(m1,m2);

    m1.setConstants(2);
    set(m2, {2,2,2,2});
    EXPECT_EQ(m1, m2);

    m1.setIdentity(3, 3);
    m2.resize(3,3);
    set(m2, {1,0,1,0});
    EXPECT_EQ(m1, m2);

    m1.setZero(4,4);
    m2.resize(4,4);
    set(m2, {0,0,0,0});
    EXPECT_EQ(m1,m2);

    m1.setConstants(5, 5, 2);
    m2.resize(5,5);
    set(m2, {2,2,2,2});
    EXPECT_EQ(m1, m2);
}

// Тест операции инверсии матрицы
TEST(MatrixTest, Inverse) {
    Matrix a(2, 2);
    set(a, {4, 7, 2, 6});
    Matrix inv = a.inverse();
    Matrix identity = a * inv;
    EXPECT_NEAR(identity.coeffRef(0, 0), 1.0, EPS);
    EXPECT_NEAR(identity.coeffRef(0, 1), 0.0, EPS);
    EXPECT_NEAR(identity.coeffRef(1, 0), 0.0, EPS);
    EXPECT_NEAR(identity.coeffRef(1, 1), 1.0, EPS);
}

// Тест транспонирования
TEST(MatrixTest, Transpose) {
    Matrix a(2, 3);
    set(a, {1, 2, 3, 4, 5, 6});
    Matrix at = a.transpose();
    EXPECT_EQ(at.rows(), 3);
    EXPECT_EQ(at.cols(), 2);
    EXPECT_DOUBLE_EQ(at.coeffRef(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(at.coeffRef(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(at.coeffRef(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(at.coeffRef(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(at.coeffRef(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(at.coeffRef(2, 1), 6.0);
}

// Тест вычисления детерминанта
TEST(MatrixTest, Determinant) {
    Matrix a(2, 2);
    set(a,{4, 3, 3, 2});
    EXPECT_DOUBLE_EQ(a.det(), -1.0);
}

// TEST Invalids

TEST(MatrixTest, INVALIDS) {
    Matrix a(2, 2);
    Matrix b(3, 3);
    Matrix c(3, 3);

    c = a*b;
    EXPECT_DOUBLE_EQ(c.isValid(), false);
}