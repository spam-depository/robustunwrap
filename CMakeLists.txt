cmake_minimum_required(VERSION 3.17)
project(RobustUnwrap
        VERSION 2023.0
        DESCRIPTION "Robust MRI phase unwrapping"
        LANGUAGES CXX C
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED yes)

add_subdirectory(src)
add_subdirectory(interfaces)
