cmake_minimum_required(VERSION 3.17)

option(PYTHON_INT "Enable building Python interface" OFF)
option(MATLAB_INT "Enable building Matlab interface" OFF)

if(PYTHON_INT)
    add_subdirectory(python)
endif (PYTHON_INT)

if(MATLAB_INT)
    add_subdirectory(matlab)
endif (MATLAB_INT)