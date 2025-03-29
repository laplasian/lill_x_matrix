//
// Created by user on 27.03.2025.
//

#include <gtest/gtest.h>
#include "matrix.h"  // Подключаем заголовочный файл с реализацией класса Matrix

// Тест конструктора и базовых операций
TEST(MatrixTest, ConstructorAndBasicOperations) {
    Matrix m(2, 2);
    EXPECT_EQ(m.rows(), 2);
    EXPECT_EQ(m.cols(), 2);

    m(0, 0) = 1.0;
    m(0, 1) = 2.0;
    m(1, 0) = 3.0;
    m(1, 1) = 4.0;
    EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(m(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(m(1, 1), 4.0);
}

// Тест конструктора копирования
TEST(MatrixTest, CopyConstructor) {
    Matrix a(2, 2, {1, 2, 3, 4});
    Matrix b(a);
    EXPECT_EQ(a, b);
}

// Тест оператора присваивания
TEST(MatrixTest, AssignmentOperator) {
    Matrix a(2, 2, {1, 2, 3, 4});
    Matrix b = a;
    EXPECT_EQ(a, b);
}

// Тест операции сложения
TEST(MatrixTest, Addition) {
    Matrix a(2, 2, {1, 2, 3, 4});
    Matrix b(2, 2, {4, 3, 2, 1});
    Matrix c = a + b;
    EXPECT_DOUBLE_EQ(c(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(c(0, 1), 5.0);
    EXPECT_DOUBLE_EQ(c(1, 0), 5.0);
    EXPECT_DOUBLE_EQ(c(1, 1), 5.0);
}

// Тест операции умножения
TEST(MatrixTest, Multiplication) {
    Matrix a(2, 2, {1, 2, 3, 4});
    Matrix b(2, 2, {2, 0, 1, 2});
    Matrix c = a * b;
    EXPECT_DOUBLE_EQ(c(0, 0), 4.0);
    EXPECT_DOUBLE_EQ(c(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(c(1, 0), 10.0);
    EXPECT_DOUBLE_EQ(c(1, 1), 8.0);
}

// Тест операции инверсии матрицы
TEST(MatrixTest, Inverse) {
    Matrix a(2, 2, {4, 7, 2, 6});
    Matrix inv = a.inverse();
    Matrix identity = a * inv;
    EXPECT_NEAR(identity(0, 0), 1.0, EPS);
    EXPECT_NEAR(identity(0, 1), 0.0, EPS);
    EXPECT_NEAR(identity(1, 0), 0.0, EPS);
    EXPECT_NEAR(identity(1, 1), 1.0, EPS);
}

// Тест транспонирования
TEST(MatrixTest, Transpose) {
    Matrix a(2, 3, {1, 2, 3, 4, 5, 6});
    Matrix at = a.transpose();
    EXPECT_EQ(at.rows(), 3);
    EXPECT_EQ(at.cols(), 2);
    EXPECT_DOUBLE_EQ(at(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(at(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(at(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(at(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(at(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(at(2, 1), 6.0);
}

// Тест вычисления детерминанта
TEST(MatrixTest, Determinant) {
    Matrix a(2, 2, {4, 3, 3, 2});
    EXPECT_DOUBLE_EQ(a.det(), -1.0);
}
