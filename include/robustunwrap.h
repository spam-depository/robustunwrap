#pragma once

#include <cstddef>

void unwrap_helper(ptrdiff_t seedx,
                   ptrdiff_t seedy,
                   ptrdiff_t seedz,
                   ptrdiff_t UNWRAPBINS,
                   const ptrdiff_t *dim_,
                   ptrdiff_t sze_,
                   ptrdiff_t m_bsx_,
                   ptrdiff_t m_bsy_,
                   ptrdiff_t m_bsz_,
                   const double *phase_,
                   const double *mag_,
                   double *unwrapped_);