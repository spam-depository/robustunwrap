#pragma once

#include <cstddef>

typedef ptrdiff_t index_t; // for specifying indices and offsets into the images

void unwrap_helper(index_t seedx,
                   index_t seedy,
                   index_t seedz,
                   index_t num_unwrapbins,
                   const index_t *dim_,
                   index_t sze_,
                   index_t m_bsx_,
                   index_t m_bsy_,
                   index_t m_bsz_,
                   const double *phase_,
                   const double *mag_,
                   double *unwrapped_);