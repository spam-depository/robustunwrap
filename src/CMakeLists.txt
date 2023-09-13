file(GLOB_RECURSE
        HEADER_LIST
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        "${Subsample_SOURCE_DIR}/include/*.h*")

add_library(robustunwrap_library
        STATIC
        robustunwrap.cpp
)

target_include_directories(robustunwrap_library PUBLIC "${RobustUnwrap_SOURCE_DIR}/include")

target_compile_features(robustunwrap_library PUBLIC cxx_std_17)
target_compile_options(robustunwrap_library PRIVATE -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable)
set_target_properties(robustunwrap_library
        PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        INTERFACE_POSITION_INDEPENDENT_CODE ON
        CXX_VISIBILITY_PRESET hidden
        C_VISIBILITY_PRESET hidden
        CUDA_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
)

source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${HEADER_LIST})

install (TARGETS robustunwrap_library
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/include" COMPONENT Development
)
