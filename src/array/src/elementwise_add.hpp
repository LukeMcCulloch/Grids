#pragma once


#include "../common/common.h"
#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <time.h>
#include <cuda_runtime.h>

//#include "etops1.hpp"

//template <typename T, typename OP1, typename OP2>
void  A_Add::sum1(int nrows, int ncols) const{

    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            val_(i,j) = op1(i,j) + op2(i,j);
        }
    }

    }