cmake_minimum_required(VERSION 3.10)

project(fmi_svd CXX)

set(CMAKE_CXX_STANDARD 17)

ADD_DEFINITIONS(-O2)

add_subdirectory(./fmi)

include_directories(./inc ${FMI_INCLUDE_DIRS} /opt/homebrew/Cellar/openblas/0.3.24/include /opt/homebrew/include)

add_library(fmipca src/fmi_wrapper.cpp src/kiss.cpp src/mmio_dense.cpp
            src/mmio.cpp src/svd_utils.cpp src/utils.cpp src/svd_algs_fmi.cpp)

link_libraries(fmipca)

add_executable(fmi_svd src/fmi_svd.cpp)
add_executable(fmi_pca src/fmi_pca.cpp)
add_executable(fmi_svd_bench src/fmi_svd_bench.cpp)

target_link_libraries(fmi_svd PRIVATE FMI "-L/opt/homebrew/Cellar/lapack/3.12.0/lib -L/opt/homebrew/Cellar/openblas/0.3.24/lib -llapacke -lopenblas")
target_link_libraries(fmi_pca PRIVATE FMI "-L/opt/homebrew/Cellar/lapack/3.12.0/lib -L/opt/homebrew/Cellar/openblas/0.3.24/lib -llapacke -lopenblas")
target_link_libraries(fmi_svd_bench PRIVATE FMI "-L/opt/homebrew/Cellar/lapack/3.12.0/lib -L/opt/homebrew/Cellar/openblas/0.3.24/lib -llapacke -lopenblas")

target_include_directories(fmi_svd PRIVATE ${FMI_INCLUDE_DIRS})
target_include_directories(fmi_pca PRIVATE ${FMI_INCLUDE_DIRS})
target_include_directories(fmi_svd_bench PRIVATE ${FMI_INCLUDE_DIRS})
