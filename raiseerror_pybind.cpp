#include <stdexcept>
#include "raiseerror.h"

void raiseerror_api(char const *msg) {
    throw std::runtime_error(msg);
};