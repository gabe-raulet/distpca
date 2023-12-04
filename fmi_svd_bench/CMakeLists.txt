cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)
project(fmi_svd_bench LANGUAGES CXX)
add_subdirectory(./fmi)
add_definitions(-O2)

find_package(aws-lambda-runtime REQUIRED)
find_package(AWSSDK COMPONENTS s3)
#find_package(LAPACKE REQUIRED)

include(FindBLAS)


include_directories(./inc ${FMI_INCLUDE_DIRS})
#include_directories(./inc ${FMI_INCLUDE_DIRS} /opt/lapack-3.12.0/CBLAS/include /opt/lapack-3.12.0/LAPACKE/include)

add_library(fmipca src/fmi_wrapper.cpp src/kiss.cpp src/mmio_dense.cpp
            src/mmio.cpp src/svd_utils.cpp src/utils.cpp src/svd_algs_fmi.cpp)

link_libraries(fmipca)

add_executable(fmi_svd_bench fmi_svd_bench.cpp)
target_link_libraries(fmi_svd_bench PUBLIC FMI AWS::aws-lambda-runtime ${AWSSDK_LINK_LIBRARIES})
target_link_libraries(fmi_svd_bench PRIVATE lapacke "-lcblas")
#target_link_libraries(fmi_svd_bench PUBLIC lapacke lapack cblas blas m)
#target_link_libraries(fmi_svd_bench PUBLIC "/opt/lapack-3.12.0/liblapack.a")
#target_link_libraries(fmi_svd_bench PUBLIC "/opt/lapack-3.12.0/libcblas.a")
#target_link_libraries(fmi_svd_bench PUBLIC "/opt/lapack-3.12.0/liblapacke.a")
#target_link_libraries(fmi_svd_bench PUBLIC "-llapacke -llapack -lblas")
#target_link_libraries(fmi_svd_bench PUBLIC FMI AWS::aws-lambda-runtime "/opt/lapack-3.12.0/liblapacke.a /opt/lapack-3.12.0/libcblas.a -L/opt/lapack-3.12.0/ -llapacke -lcblas")
target_include_directories(fmi_svd_bench PRIVATE ${FMI_INCLUDE_DIRS})
target_compile_features(fmi_svd_bench PRIVATE "cxx_std_11")
target_compile_options(fmi_svd_bench PRIVATE "-O2")

aws_lambda_package_target(fmi_svd_bench)
