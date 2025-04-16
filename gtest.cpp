#include <cmath>
#include <gtest/gtest.h>
#include "Matrix.hpp"

Matrix construct(std::vector<std::vector<double>> data) {
    size_t raws = data.size();
    size_t cols = data[0].size();
    for (const auto & i : data) {
        if (i.size() != cols) {
            throw std::out_of_range("invalid matrix data");
        }
    }
    Matrix mat(raws, cols);

    for (int i = 0; i < raws; i++) {
        for (int j = 0; j < cols; j++) {
            mat.coeffRef(i, j) = data[i][j];
        }
    }

    return mat;
}

bool operator==(const Matrix &a, const Matrix &b) {
    if (!a.isValid() || !b.isValid() || a.cols()!= b.cols() || a.rows() != b.rows()) return false;
    for (int i = 0; i < a.cols()*a.rows(); i++) {
        if (a.data()[i] != b.data()[i]) return false;
    }
    return true;
}

// === Constructors ===
TEST(MatrixTest, DefaultConstructor) {
    Matrix m;
    EXPECT_FALSE(m.isValid());
}

TEST(MatrixTest, ConstructorWithRowsCols) {
    Matrix m(2, 3);
    EXPECT_EQ(m.rows(), 2);
    EXPECT_EQ(m.cols(), 3);
    EXPECT_TRUE(m.isValid());
}

TEST(MatrixTest, CopyConstructorValid) {
    auto m1 = construct({{1, 2}, {3, 4}});
    Matrix m2(m1);
    EXPECT_TRUE(m2.isValid());
    EXPECT_EQ(m1, m2);
}

TEST(MatrixTest, AssignmentOperator) {
    auto m1 = construct({{1, 2}, {3, 4}});
    Matrix m2;
    m2 = m1;
    EXPECT_EQ(m1, m2);
}

// === Operators ===
TEST(MatrixTest, OperatorAdditionValid) {
    auto a = construct({{1, 2}, {3, 4}});
    auto b = construct({{5, 6}, {7, 8}});
    auto c = a + b;
    EXPECT_EQ(c, construct({{6, 8}, {10, 12}}));
}

TEST(MatrixTest, OperatorAdditionInvalid) {
    auto a = construct({{1, 2}, {3, 4}});
    auto b = construct({{1, 2}});
    auto c = a + b;
    EXPECT_FALSE(c.isValid());
}

TEST(MatrixTest, OperatorSubtractionValid) {
    auto a = construct({{5, 6}, {7, 8}});
    auto b = construct({{1, 2}, {3, 4}});
    auto c = a - b;
    EXPECT_EQ(c, construct({{4, 4}, {4, 4}}));
}

TEST(MatrixTest, OperatorMultiplicationValid) {
    auto a = construct({{1, 2}, {3, 4}});
    auto b = construct({{2, 0}, {1, 2}});
    auto c = a * b;
    EXPECT_EQ(c, construct({{4, 4}, {10, 8}}));
}

TEST(MatrixTest, OperatorScalarMultiplication) {
    auto a = construct({{1, 2}, {3, 4}});
    auto c = a * 2.0;
    EXPECT_EQ(c, construct({{2, 4}, {6, 8}}));
}

TEST(MatrixTest, OperatorScalarDivisionValid) {
    auto a = construct({{2, 4}, {6, 8}});
    auto c = a / 2.0;
    EXPECT_EQ(c, construct({{1, 2}, {3, 4}}));
}

TEST(MatrixTest, OperatorScalarDivisionByZero) {
    auto a = construct({{2, 4}, {6, 8}});
    EXPECT_THROW(a / 0.0, std::runtime_error);
}

TEST(MatrixTest, OperatorMulAssignValid) {
    auto a = construct({{1, 2}, {3, 4}});
    auto b = construct({{2, 0}, {1, 2}});
    a *= b;
    EXPECT_EQ(a, construct({{4, 4}, {10, 8}}));
}

TEST(MatrixTest, OperatorAddAssignValid) {
    auto a = construct({{1, 1}, {1, 1}});
    auto b = construct({{2, 2}, {2, 2}});
    a += b;
    EXPECT_EQ(a, construct({{3, 3}, {3, 3}}));
}

TEST(MatrixTest, OperatorSubAssignValid) {
    auto a = construct({{3, 3}, {3, 3}});
    auto b = construct({{1, 1}, {1, 1}});
    a -= b;
    EXPECT_EQ(a, construct({{2, 2}, {2, 2}}));
}

TEST(MatrixTest, OperatorMulAssignScalar) {
    auto a = construct({{1, 2}, {3, 4}});
    a *= 3.0;
    EXPECT_EQ(a, construct({{3, 6}, {9, 12}}));
}

TEST(MatrixTest, OperatorDivAssignScalarValid) {
    auto a = construct({{4, 6}, {8, 10}});
    a /= 2.0;
    EXPECT_EQ(a, construct({{2, 3}, {4, 5}}));
}

// === Transpose ===
TEST(MatrixTest, Transpose) {
    auto a = construct({{1, 2, 3}, {4, 5, 6}});
    auto t = a.transpose();
    EXPECT_EQ(t, construct({{1, 4}, {2, 5}, {3, 6}}));
}

// === Determinant ===
TEST(MatrixTest, Determinant2x2) {
    auto a = construct({{1, 2}, {3, 4}});
    EXPECT_NEAR(a.det(), -2.0, 1e-6);
}

TEST(MatrixTest, DeterminantNonSquare) {
    auto a = construct({{1, 2, 3}, {4, 5, 6}});
    EXPECT_TRUE(std::isnan(a.det()));
}

// === Inverse ===
TEST(MatrixTest, Inverse2x2) {
    auto a = construct({{4, 7}, {2, 6}});
    auto inv = a.inverse();
    auto identity = a * inv;
    for (int i = 0; i < identity.rows(); i++) {
        for (int j = 0; j < identity.cols(); j++) {
            if (i == j) EXPECT_NEAR(identity.coeffRef(i, j), 1.0, 1e-6);
            else EXPECT_NEAR(identity.coeffRef(i, j), 0.0, 1e-6);
        }
    }
}

TEST(MatrixTest, InverseSingular) {
    auto a = construct({{1, 2}, {2, 4}});
    auto inv = a.inverse();
    EXPECT_FALSE(inv.isValid());
}

// === Set Identity/Zero/Constants ===
TEST(MatrixTest, SetIdentity) {
    Matrix a(3, 3);
    a.setIdentity();
    EXPECT_EQ(a, construct({{1,0,0},{0,1,0},{0,0,1}}));
}

TEST(MatrixTest, SetZero) {
    Matrix a(2, 2);
    a.setZero();
    EXPECT_EQ(a, construct({{0, 0}, {0, 0}}));
}

TEST(MatrixTest, SetConstants) {
    Matrix a(2, 2);
    a.setConstants(3.14);
    EXPECT_EQ(a, construct({{3.14, 3.14}, {3.14, 3.14}}));
}

// === Resize ===
TEST(MatrixTest, ResizeValid) {
    Matrix a(2, 2);
    a.resize(3, 3);
    EXPECT_EQ(a.rows(), 3);
    EXPECT_EQ(a.cols(), 3);
}

// === Exceptions ===
TEST(MatrixTest, CoeffRefOutOfBounds) {
    Matrix a(2, 2);
    EXPECT_THROW(a.coeffRef(5, 5), std::out_of_range);
}

TEST(MatrixTest, InvalidMatrixCoeffRef) {
    Matrix a;
    EXPECT_THROW(a.coeffRef(0, 0), std::out_of_range);
}

// === Static Constructors ===
TEST(MatrixTest, StaticIdentity) {
    auto id = Matrix::identity(2, 2);
    EXPECT_EQ(id, construct({{1, 0}, {0, 1}}));
}

TEST(MatrixTest, StaticZeros) {
    auto z = Matrix::zeros(2, 2);
    EXPECT_EQ(z, construct({{0, 0}, {0, 0}}));
}

TEST(MatrixTest, StaticConstants) {
    auto c = Matrix::constants(2, 2, 42.0);
    EXPECT_EQ(c, construct({{42.0, 42.0}, {42.0, 42.0}}));
}
