#pragma once


#include <stddef.h>
#include <cassert>

// include helper class traits template to select whether to refer to an
// ``expression template node'' either ``by value'' or ``by reference.''
#include "etops1a.hpp"
#include "etmatrix.hpp"

//#include "elementwise_add.hpp"

using namespace std;

// class for objects that represent the addition of two operands
template <typename T, typename OP1, typename OP2>
class A_Add {
  private:
    typename A_Traits<OP1>::ExprRef op1;    // first operand
    typename A_Traits<OP2>::ExprRef op2;    // second operand

  public: 
    bool iscached = 0;
    bool evalme = 0;

    int storage_size, nrows, ncols;

    mutable Array<T> val_;

    // constructor initializes references to operands
    A_Add (OP1 const& a, OP2 const& b)
     : op1(a), op2(b) {
       
      //fixme: assert op1.size() == op2.size();
      //fixme: assert op1.getnrows() == op2.getnrows();
      nrows = op1.getnrows();
      ncols = op1.getncols();
      storage_size = op1.size();
      // testing only!:
      //cache();
      //iscached =1;
    }

    //void sum1(int,int) const;
    void  sum1(int nrows, int ncols) const{

      for (size_t i = 0; i < nrows; i++) {
          for (size_t j = 0; j < ncols; j++) {
            val_(i,j) = op1(i,j) + op2(i,j);
          }
      }
    }

    void  sum2(int nrows, int ncols) const{
      const int N = op1.size();
      for (size_t i = 0; i < N; i++) {
            val_(i) = op1(i) + op2(i);
          }
    }
    

    void cache() const {
      int nrows = op1.getnrows();
      int ncols = op1.getncols();
      val_.SetArray(op1.getnrows(), op2.getncols());
      op1.cache(); op2.cache();

      val_ = 0.;

      sum2(nrows,ncols);
      //sum_gpu_setup(); //helper 

      

      //val_ = op1.val() + op2.val();
      
      //iscached =1;
    }

    Array<T>& val() const{
      //return *val_;
      return *val_.expr_rep.val();
    }


    //compute sum when value requested
    T operator()(size_t i, size_t j) const{
      if (iscached==1){
        return val_(i,j);
      }
      else{
        assert(i<op1.getnrows());
        assert(j<op1.getncols());
        return op1(i,j) + op2(i,j);
      }
    }

    // size is maximum size
    size_t size() const {
        assert (op1.size()==0 || op2.size()==0
                || op1.size()==op2.size());
        return op1.size()!=0 ? op1.size() : op2.size();
    }

    
      size_t getnrows() const {
          assert (op1.getnrows()==0 || op2.getnrows()==0
                || op1.getnrows()==op2.getnrows());
          return op1.getnrows()!=0 ? op1.getnrows() : op2.getnrows();
      }
      
      size_t getncols() const {
          assert (op1.getncols()==0 || op2.getncols()==0
                || op1.getncols()==op2.getncols());
          return op1.getncols()!=0 ? op1.getncols() : op2.getncols();
      }

      // size_t getnrows() const {
      //     return op1.getnrows();
      // }
      
      // size_t getncols() const {
      //     return op1.getncols();
      // }

};

// template <typename T, typename OP1, typename OP2>
// void  A_Add::sum1(int nrows, int ncols) const{

//   for (size_t i = 0; i < nrows; i++) {
//       for (size_t j = 0; j < ncols; j++) {
//           val_(i,j) = op1(i,j) + op2(i,j);
//       }
//   }

// }




// class for objects that represent the multiplication of two operands
template <typename T, typename OP1, typename OP2>
class A_Mult {
  private:
    typename A_Traits<OP1>::ExprRef op1;    // first operand
    typename A_Traits<OP2>::ExprRef op2;    // second operand

  public:
    mutable bool iscached = 0;
    bool evalme = 0;

    mutable Array<T> val_;

    // constructor initializes references to operands
    A_Mult (OP1 const& a, OP2 const& b)
     : op1(a), op2(b) {


      // testing only!:
      //cache();
      //iscached =1;
    }


    void cache() const {
      int nrows = op1.getnrows();
      int ncols = op1.getncols();

      val_.SetArray(op1.getnrows(), op2.getncols());
      op1.cache(); op2.cache();

      val_ = 0.;

      for (size_t i = 0; i < nrows; i++) {
          for (size_t j = 0; j < ncols; j++) {
            val_(i,j) = op1(i,j) + op2(i,j);
          }
      }
      //iscached =1;
      //val_ = op1.val() * op2.val();
    }

    Array<T>& val() const {
      //return *val_;
      return *val_.expr_rep.val();
    }

    // compute product when value requested
    T operator() (size_t i, size_t j) const {
      if (iscached==1){
        return val_(i,j);
      }
      else{
        return op1(i,j) * op2(i,j);
      }
    }
    T operator() (size_t i) const {
      if (iscached==1){
        return val_(i);
      }
      else{
        return op1(i) * op2(i);
      }
    }

    // size is maximum size
    size_t size() const {
        assert (op1.size()==0 || op2.size()==0
                || op1.size()==op2.size());
        return op1.size()!=0 ? op1.size() : op2.size();
    }

    
      size_t getnrows() const {
          assert (op1.getnrows()==0 || op2.getnrows()==0
                || op1.getnrows()==op2.getnrows());
          return op1.getnrows()!=0 ? op1.getnrows() : op2.getnrows();
      }

      size_t getncols() const {
          assert (op1.getncols()==0 || op2.getncols()==0
                || op1.getncols()==op2.getncols());
          return op1.getncols()!=0 ? op1.getncols() : op2.getncols();
      }
};




// class for objects that represent the matrix multiplication of two operands
template <typename T, typename OP1, typename OP2>
class MatMult {
  private:
    typename A_Traits<OP1>::ExprRef op1;    // first operand
    typename A_Traits<OP2>::ExprRef op2;    // second operand
    


  public:
    bool iscached = 0;
    bool evalme = 1;

    mutable Array<T> val_;


    // constructor initializes references to operands
    MatMult (OP1 const& a, OP2 const& b)
    : op1(a), op2(b) {

      cache();
      iscached =1;
    }

    

    void cache() const{
      int nrows = op1.getnrows();
      int ncols = op1.getncols();
      val_.SetArray(op1.getnrows(), op2.getncols());
      op1.cache(); op2.cache();

      val_ = 0.;
      
      for (size_t i = 0; i < nrows; i++) {
          for (size_t k = 0; k < ncols; k++) {
              for (size_t j = 0; j < ncols; j++) {
                  val_(i, j) += op1(i, k) * op2(k, j);
              }
          }
      }


      //val_.set_array(matmul(op1.val() , op2.val() )) ;

      //iscached =1;

      //val_.set_array(matmul(op1,op2));
    }

    Array<T>& val()const{
      //return *val_;
      return *val_.expr_rep.val();
    }

    // compute product when value requested
    // lazy, expression template "naive" usage pattern:
    // T operator() (size_t i, size_t j) const {
    //   T temp = 0.;
    //   for (size_t k = 0; k < op1.getncols(); k++) {
    //     temp += op1(i, k) * op2(k, j);
    //   }
    //   return temp;
    // }


    //compute sum when value requested
    // T operator()(size_t i, size_t j) const{
    //   if (iscached==1){
    //     return val_(i,j);
    //   }
    //   else{
    //     T temp = 0.;
    //     for (size_t k = 0; k < op1.getncols(); k++) {
    //       temp += op1(i, k) * op2(k, j);
    //     }
    //     return temp;
    //   }
    // }

    // result is cached.  return when requested
    T operator() (size_t i, size_t j) const {
      return val_(i,j);
    }



    // size is maximum size
    size_t size() const {
        assert (op1.size()==0 || op2.size()==0
                || op1.size()==op2.size());
        return op1.size()!=0 ? op1.size() : op2.size();
    }

    
      size_t getnrows() const {
          assert (op1.getnrows()==0 || op2.getnrows()==0
                || op1.getnrows()==op2.getnrows());
          return op1.getnrows()!=0 ? op1.getnrows() : op2.getnrows();
      }

      size_t getncols() const {
          assert (op1.getncols()==0 || op2.getncols()==0
                || op1.getncols()==op2.getncols());
          return op1.getncols()!=0 ? op1.getncols() : op2.getncols();
      }
};



// class for objects that represent the matrix multiplication of two operands
template <typename T, typename OP1, typename OP2>
class MatDotVec {
  private:
    typename A_Traits<OP1>::ExprRef op1;    // first operand
    typename A_Traits<OP2>::ExprRef op2;    // second operand

    //template <typename T, Rep = Array2D<T> >
  public:
    bool iscached = 0;
    bool evalme = 1;
  
    mutable Array<T> val_;

    Array<T> result2;

    //ExprRef cache_result;
    //implement_matmul;
    // constructor initializes references to operands
    MatDotVec (OP1 const& a, OP2 const& b)
    : op1(a), op2(b) {
      //implement_matmul = matmul;
       //cache_result = matmul(op1,op2);
    }


    // auto operator=(){
    //   auto c(i,k) += op1.array(i,j) * op2.array(j,k);
    //   return c
    // }

    // compute product when value requested
    // auto operator=() const {
    //     auto result = matmul(op1 , op2) ;
    //     return result(i,j);
    // }

    // compute product when value requested
    // OP1 matmul () {
    //     return matmul(op1,op2) ;
    // }
    

    // compute product when value requested
    T operator() (size_t i) const {
      T temp = 0.;
      for (size_t k = 0; k < op1.getncols(); k++) {
        temp += op1(i, k) * op2(k);
      }
      return temp;
    }

    // template <typename C>
    // void assign_to(C&& c) const {

    //     auto& op1 = this->op1();
    //     auto& op2 = this->op2();

    // }

    // size is maximum size
    size_t size() const {
        assert (op1.size()==0 || op2.size()==0
                || op1.size()==op2.size());
        return op1.size()!=0 ? op1.size() : op2.size();
    }

    
      size_t getnrows() const {
          assert (op1.getnrows()==0 || op2.getnrows()==0
                || op1.getnrows()==op2.getnrows());
          return op1.getnrows()!=0 ? op1.getnrows() : op2.getnrows();
      }

      size_t getncols() const {
          assert (op1.getncols()==0 || op2.getncols()==0
                || op1.getncols()==op2.getncols());
          return op1.getncols()!=0 ? op1.getncols() : op2.getncols();
      }
};
