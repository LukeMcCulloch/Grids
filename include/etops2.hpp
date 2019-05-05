

#include "arrayops.hpp"



using namespace std;



// addition of two Arrays
template <typename T, typename R1, typename R2>
Array<T,A_Add<T,R1,R2> >
operator+ (Array<T,R1> const& a, Array<T,R2> const& b) {
    //cout << "adding two arrays in etops2.hpp\n" ;
    return Array<T,A_Add<T,R1,R2> >
           (A_Add<T,R1,R2>(a.rep(),b.rep()));
}


// multiplication of two Arrays
template <typename T, typename R1, typename R2>
Array<T, A_Mult<T,R1,R2> >
operator* (Array<T,R1> const& a, Array<T,R2> const& b) {
    //cout << "multiplying two arrays in etops2.hpp\n" ;
    return Array<T,A_Mult<T,R1,R2> >
           (A_Mult<T,R1,R2>(a.rep(), b.rep()));
}

// multiplication of scalar and Array
template <typename T, typename R2>
Array<T, A_Mult<T,A_Scalar<T>,R2> >
operator* (T const& s, Array<T,R2> const& b) {
    //cout << "multiplying an array and a scalar in etops2.hpp\n" ;
    return Array<T,A_Mult<T,A_Scalar<T>,R2> >
           (A_Mult<T,A_Scalar<T>,R2>(A_Scalar<T>(s), b.rep()));
}

// multiplication of Array and scalar
// addition of scalar and Array
// addition of Array and scalar
//...


// // matmul two Arrays
template <typename T, typename R1, typename R2>
Array<T, MatMult<T,R1,R2> >
matmul (Array<T,R1> const& a, Array<T,R2> const& b) {
    //cout << "multiplying two arrays in etops2.hpp\n" ;
    return Array<T,MatMult<T,R1,R2> >
           (MatMult<T,R1,R2>(a.rep(), b.rep()));
}

// multiplication of two Arrays
// template <typename T, typename R1, typename R2>
// Array<T, MatMult<T,R1,R2> >
// operator* (Array<T,R1> const& a, Array<T,R2> const& b) {
//     //cout << "multiplying two arrays in etops2.hpp\n" ;
//     return Array<T,MatMult<T,R1,R2> >
//            (MatMult<T,R1,R2>(a.rep(), b.rep()));
// }