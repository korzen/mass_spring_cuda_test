cmake_minimum_required(VERSION 2.8)

project(fast_mass_spring)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#set OpenMp useable
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11 -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED -Wno-deprecated")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#message("xxxx"+${CMAKE_C_FLAGS})

#make eigen available
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3/") 
include_directories( "${EIGEN3_INCLUDE_DIR}" )

#make spectra available
set(SPECTRA_INCLUDE_DIR "/usr/local/include/spectra-master/include/")
include_directories( "${SPECTRA_INCLUDE_DIR}" )

#make boost available
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.58.0 COMPONENTS date_time filesystem system program_options REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message(${Boost_LIBRARIES})
else()
  message("--[info] can not find Boost")
endif()

#make eigen SuiteSparse available
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR};${CAMKE_MODULE_PATH}")
find_package(SuiteSparse)
include_directories(${SuiteSparse_INCLUDE_DIRS})

#get file for linking 
aux_source_directory(. src)
file(GLOB header "*.h")
file(GLOB cuda_file "*.cu")
message("src::" ${src})
message("header file::" ${header})
message("cuda file::" ${cuda_file})

#compiler target
#link_directories("/usr/lib/x86_64-linux-gnu/")

#find_package(CUDA)
FIND_PACKAGE(CUDA ${CUDA_VERSION} REQUIRED)
STRING (TOLOWER ${CMAKE_SYSTEM_NAME} LOWER_SYSTEM_NAME)
SET(CUDA_TARGET_INCLUDE  ${CUDA_TOOLKIT_ROOT_DIR}-${CUDA_VERSION}/targets/${CMAKE_HOST_SYSTEM_PROCESSOR}-${LOWER_SYSTEM_NAME}/include)
message(${CUDA_INCLUDE_DIRS})

SET(CUDA_NVCC_FLAGS --compiler-options "-std=c++11")

#generate exe
cuda_add_executable(fast_mass_spring ${src} ${header} ${cuda_file} )
target_link_libraries(fast_mass_spring GLU glut GL ${Boost_LIBRARIES} ${SuiteSparse_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} )
