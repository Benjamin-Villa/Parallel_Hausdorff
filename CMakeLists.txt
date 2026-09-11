cmake_minimum_required(VERSION 3.22)

# Find SYCL compiler (icpx or clang++)
find_program(ICPX_COMPILER icpx)
if(ICPX_COMPILER)
    set(CMAKE_CXX_COMPILER ${ICPX_COMPILER})
else()
    set(CMAKE_CXX_COMPILER "clang++")
endif()

project(Parallel_HSD_K2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configure CPU target based on compiler (spir64_x86_64 for icpx, native_cpu for clang++)
if(CMAKE_CXX_COMPILER MATCHES "icpx")
    set(INTEL_CPU_TARGET "spir64_x86_64")
else()
    set(INTEL_CPU_TARGET "native_cpu")
endif()

set(CPU_SYCL_FLAGS
    -fsycl
    -fsycl-targets=${INTEL_CPU_TARGET}
    -O3
)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/seq
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/K2tree
    ${CMAKE_CURRENT_SOURCE_DIR}/Util
    ${CMAKE_CURRENT_SOURCE_DIR}/libcds/include
)

# Compile libcds static library
file(GLOB_RECURSE LIBCDS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/libcds/src/*.cpp
)
add_library(cds_par STATIC ${LIBCDS_SRC})
target_include_directories(cds_par PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libcds/include)

# Executable with sequential baseline algorithms and parallel HDKP1
add_executable(Parallel_HSD_K2
    src/main.cpp
    src/ParallelHDKP1.cpp
    src/seq/NNMAX.cpp
    src/seq/IsCandidate.cpp
    src/seq/HDKP1.cpp
    src/seq/HDKP2.cpp
    K2tree/K2tree.cpp
    Util/Factory.cpp
    Util/Point.cpp
    Util/TimeMesure.cpp
    utils.cpp
)

target_compile_options(Parallel_HSD_K2 PRIVATE ${CPU_SYCL_FLAGS})
target_link_options(Parallel_HSD_K2 PRIVATE ${CPU_SYCL_FLAGS})
target_link_libraries(Parallel_HSD_K2 PRIVATE cds_par)
