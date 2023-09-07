#include <stddef.h>

void unwrap_helper(const ptrdiff_t seedx, const ptrdiff_t seedy, const ptrdiff_t seedz,
                   const ptrdiff_t UNWRAPBINS,
                   const ptrdiff_t *dim_,
                   const ptrdiff_t sze_,
                   const ptrdiff_t m_bsx_,
                   const ptrdiff_t m_bsy_,
                   const ptrdiff_t m_bsz_,
                   const double *phase_,
                   const double *mag_,
                   double *unwrapped_);