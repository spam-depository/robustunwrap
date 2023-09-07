#include <stddef.h>

void unwrap_helper(int seedx, int seedy, int seedz,
                   long UNWRAPBINS,
                   ptrdiff_t *dim_,
                   size_t sze_,
                   ptrdiff_t m_bsx_,
                   ptrdiff_t m_bsy_,
                   ptrdiff_t m_bsz_,
                   double *phase_,
                   double *mag_,
                   double *unwrapped_);