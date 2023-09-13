macro( build_mex name func_src private_deps)
    message(STATUS "matlab module ${name} will be installed to ${CMAKE_INSTALL_FULL_LIBDIR}")

    #build the library as -shared or as -static
    matlab_add_mex(NAME ${name} MODULE SRC "${func_src}" LINK_TO ${private_deps} R2018a)
    target_compile_features(${name} PUBLIC cxx_std_17)

    # install libraries
    # TODO: the EXPORT isis-export doesn't make sense
    install (TARGETS ${name}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
    )
endmacro()

find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER)

build_mex(robustunwrap "robustunwrap_matlab.cpp;raiseerror_matlab.cpp" robustunwrap_library)
