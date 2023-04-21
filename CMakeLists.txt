cmake_minimum_required(VERSION 3.0.0)
project(UmeyamaAlignment VERSION 0.1.0)

include(CTest)
enable_testing()

add_library(UmeyamaAlignment UmeyamaAlignment.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
