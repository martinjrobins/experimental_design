cmake_minimum_required(VERSION 2.8.3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
project(sinusoidal_voltammetry)
FIND_PACKAGE(PythonLibs)
#FIND_PACKAGE(SUNDIALS)
FIND_PACKAGE(Boost COMPONENTS python REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${SUNDIALS_INCLUDE_DIRS})
set(BUILD_SHARED_LIBS ON)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(source_files 
    src/e_surface.cpp
    )
set(header_files 
    src/sinusoidal_voltammetry.hpp
    src/utilities.hpp
    )

add_library(sinusoidal_voltammetry ${source_files} ${header_files})
target_link_libraries(sinusoidal_voltammetry ${Boost_LIBRARIES} ${SUNDIALS_LIBRARIES})

add_library(exp_design src/python_wrapper.cpp)
target_link_libraries(exp_design sinusoidal_voltammetry ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

set_target_properties(exp_design PROPERTIES PREFIX "" )

find_package(CxxTest)
if (CXXTEST_FOUND)
    enable_testing()

    include_directories(${CXXTEST_INCLUDE_DIR} ${sinusoidal_voltammetry_SOURCE_DIR})
    add_subdirectory(tests)
endif()

