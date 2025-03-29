#include <gtest/gtest.h>
#include "matrix.h"  // Подключаем заголовочный файл с реализацией класса Matrix


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}