
//=================================
// include guard
#ifndef __MATRIX_TEMPLATE_INCLUDED__
#define __MATRIX_TEMPLATE_INCLUDED__


#include <stddef.h>
#include <cassert>
#include "array_template.hpp"

template <typename T, typename Rep = Array2D<T> >
class Array {
  private:
    Rep expr_rep;   // (access to) the data of the array

  public:
    bool iscached = 1;
    bool evalme = 0;
    //Rep expr_rep; 
    // create array with initial size
    explicit Array (size_t s, size_t t)
     : expr_rep(s,t) {
    }

    // create array from possible representation
    Array (Rep const& rb)
     : expr_rep(rb) {
    }

    // default constructor:
    Array(){}



    void SetArray(size_t s, size_t t){
        expr_rep.buildWithParameters(s,t);
    }

    // MatMult is friend:
    //template < T, typename OP1, typename OP2>
    //friend class MatMult;

    
    // evaluate expressions:
    void evaluate_set(){

    }

    void set_array(const Rep& val ){
        expr_rep = val;
    }


    // assignment operator for same type
    Array& operator=(Array const& b) { 
        assert(size()==b.size());
        int nrows = getnrows();
        int ncols = getncols();
        for (int i = 0; i<nrows; ++i) {
            for (int j = 0; j<ncols; ++j) {
                expr_rep(i,j) = b(i,j);
            }
        }
        return *this;
    }

    Array& operator=(const T a) {
        int nrows = getnrows();
        int ncols = getncols();
        for (int i = 0; i<nrows; ++i) {
            for (int j = 0; j<ncols; ++j) {
                expr_rep(i,j) = a;
            }
        }
        return *this;
    }

    // assignment operator for arrays of different type
    template<typename T2, typename Rep2>
    Array& operator=(Array<T2, Rep2> const& b) { 
        assert(size()==b.size());
        int nrows = getnrows();
        int ncols = getncols();
        for (int i = 0; i<nrows; ++i) {
            for (int j = 0; j<ncols; ++j) {
                expr_rep(i,j) = b(i,j);
            }
        }
        return *this;
    }


    // size is size of represented data
    size_t getnrows() const { return expr_rep.getnrows();}
    size_t getncols() const { return expr_rep.getncols();}
    size_t size() const { return expr_rep.size();}


    void cache() const {}

    // Array& val() const {
    //     return *this;
    // }

    Rep& val() const {
        return *expr_rep.val();
    }

    // parens ops:
    T operator()(size_t i, size_t j) const {
        assert(i<getnrows());
        assert(j<getncols());
        //return expr_rep[i*getncols() + j];
        return expr_rep(i,j);
    }
    // self uses this!
    T& operator()(size_t i, size_t j) {
        assert(i<getnrows());
        assert(j<getncols());
        //return expr_rep.array[i*getncols() + j];
        return expr_rep(i,j);
    }
    
    T operator()(size_t i) const {
        assert(i<size());
        //return expr_rep[i*getncols() + j];
        return expr_rep(i);
    }
    // self uses this!
    T& operator()(size_t i) {
        assert(i<size());
        //return expr_rep.array[i*getncols() + j];
        return expr_rep(i);
    }


    // index operator for constants and variables
    // TODO:
    // index operator for constants and variables
    // T operator[] (size_t idx) const {
    //     assert(idx<size());
    //     return expr_rep[idx];
    // }
    // T& operator[] (size_t idx) {
    //     assert(idx<size());
    //     return expr_rep[idx];
    // }

    // return what the array currently represents
    Rep const& rep() const {
        return expr_rep;
    }
    Rep& rep() {
        return expr_rep;
    }
};

#endif //__MATRIX_TEMPLATE_INCLUDED__