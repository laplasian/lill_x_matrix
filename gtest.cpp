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

// Тест конструктора и базовых операций
TEST(MatrixTest, ConstructorAndBasicOperations) {
    Matrix m(2, 2);
    EXPECT_EQ(m.rows(), 2);
    EXPECT_EQ(m.cols(), 2);

    m.coeffRef(0,0) = 1.0;
    m.coeffRef(0, 1) = 2.0;
    m.coeffRef(1, 0) = 3.0;
    m.coeffRef(1, 1) = 4.0;
    EXPECT_DOUBLE_EQ(m.coeffRef(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m.coeffRef(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(m.coeffRef(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(m.coeffRef(1, 1), 4.0);
}

// Тест конструктора копирования
TEST(MatrixTest, CopyConstructor) {
    Matrix a(2, 2);
    set(a, {1,2,3,4});
    Matrix b(a);
    EXPECT_EQ(a, b);
}

// Тест оператора присваивания
TEST(MatrixTest, AssignmentOperator) {
    Matrix a(2, 2);
    set(a, {1,2,3,4});
    Matrix b = a;
    EXPECT_EQ(a, b);
}

// Тест операции сложения
TEST(MatrixTest, Addition) {
    Matrix a(2, 2);
    set(a, {1, 2, 3, 4});
    Matrix b(2, 2);
    set(b, {4, 3, 2, 1});
    Matrix c = a + b;
    EXPECT_DOUBLE_EQ(c.coeffRef(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(c.coeffRef(0, 1), 5.0);
    EXPECT_DOUBLE_EQ(c.coeffRef(1, 0), 5.0);
    EXPECT_DOUBLE_EQ(c.coeffRef(1, 1), 5.0);
}

// Тест операции умножения
TEST(MatrixTest, Multiplication) {
    Matrix a(2, 2);
    set(a, {1, 2, 3, 4});
    Matrix b(2, 2);
    set(b, {2, 0, 1, 2});
    Matrix c = a * b;
    EXPECT_DOUBLE_EQ(c.coeffRef(0, 0), 4.0);
    EXPECT_DOUBLE_EQ(c.coeffRef(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(c.coeffRef(1, 0), 10.0);
    EXPECT_DOUBLE_EQ(c.coeffRef(1, 1), 8.0);
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
