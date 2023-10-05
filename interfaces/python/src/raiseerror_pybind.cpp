#include <stdexcept>
#include "raiseerror.h"

void raiseerror(char const *msg) {
    throw std::runtime_error(msg);
};