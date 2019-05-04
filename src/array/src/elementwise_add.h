



    void  A_Add::sum1(int nrows, int ncols) const{

        for (size_t i = 0; i < nrows; i++) {
            for (size_t j = 0; j < ncols; j++) {
              val_(i,j) = op1(i,j) + op2(i,j);
            }
        }
  
      }