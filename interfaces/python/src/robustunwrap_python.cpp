/*
 This is the seperated Python Pybind interface for

 3D Unwrapping algorithm
 Rhodri Cusack 2000-2006
 Algorithm described in 
 Cusack, R. & Papadakis, N. (2002) 
    "New robust 3-D phase unwrapping algorithms: application to magnetic field mapping and undistorting echoplanar images."
    Neuroimage. 2002 Jul;16(3 Pt 1):754-64.
 Distributed under MIT License
 Comments to cusackrh@tcd.ie
 Version 3.00 Oct 2006 adapted for matlab

 Refactoring done by Podranski, K. (Sep 2023)

 Compile with:
 c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) robustunwrap_python.cpp robustunwrap.cpp raiseerror_pybind.cpp -o robustunwrap$(python3-config --extension-suffix)
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "robustunwrap.h"
#include "raiseerror.h"

namespace py = pybind11;

// TODO: Not sure if all the array/list types here could/should be a const ref
py::array_t<double, py::array::f_style | py::array::forcecast> robustunwrap(py::list seed, py::array_t<double, py::array::f_style | py::array::forcecast> phaseData, py::array_t<double, py::array::f_style | py::array::forcecast> magnitudeData, const int numunwrapbins = 10000) {
    const size_t nDims = 3;

    if (seed.size() != 3) {
        throw std::runtime_error("seed must be an index into the 3D data array consisting of 3 ints.");
    }

    const ptrdiff_t seedX = py::cast<int>(seed[0]); // maybe a bit naive. might not work that way.
    const ptrdiff_t seedY = py::cast<int>(seed[1]);
    const ptrdiff_t seedZ = py::cast<int>(seed[2]);

    
    const ptrdiff_t *dims = phaseData.shape();
    const ptrdiff_t size = phaseData.size();

    if (seedX < 0 || seedX >= dims[0] || seedY < 0 || seedY >= dims[1] || seedZ < 0 || seedZ >= dims[2]) {
        throw std::runtime_error("The seed specified was outside the matrix bounds.");
    }

    if (magnitudeData.ndim() != nDims || phaseData.ndim() != nDims) {
        throw std::runtime_error("Number of dimensions must be three");
    } else if (!std::equal(dims, dims + nDims, phaseData.shape())) {
        throw std::runtime_error("magnitudeData and phaseData shapes must match");
    }
    /* The input must be a noncomplex scalar double.*/
    // else if (magnitudeData.dtype != ??? || phaseData.dtype != ???) { ...

    /* define strides */
    ptrdiff_t m_bsx = 1;
    ptrdiff_t m_bsy = dims[0];
    ptrdiff_t m_bsz = dims[0] * dims[1];

    // Negate input as low polefield values unwrapped first
    const double *magnitudeInput = magnitudeData.data();
    auto *magnitude = new double[size];
    for (long i = 0; i < size; i++)
        magnitude[i] = -magnitudeInput[i];

    // Create matrix for the return argument.
    auto unwrappedArray = py::array_t<double, py::array::f_style | py::array::forcecast>({dims[0], dims[1], dims[2]}); // not very elegant
    double *unwrapped = static_cast<double *>(unwrappedArray.request().ptr);

    unwrap_helper(seedX, seedY, seedZ, numunwrapbins, dims, size, m_bsx, m_bsy, m_bsz,
                  phaseData.data(), magnitude, unwrapped);

    delete[] magnitude;

    return unwrappedArray;
}

PYBIND11_MODULE(robustunwrap_python_library, m) {
    m.doc() = "attempt to interface robustunwrap (Cusack & Papadakis 2002) with pybind";
    m.def("robustunwrap", &robustunwrap, "robustunwrap (Cusack & Papadakis 2002)",
        py::arg("seed"), py::arg("phaseData"), py::arg("magnitudeData"),
        py::arg("numunwrapbins") = 10000);
}
