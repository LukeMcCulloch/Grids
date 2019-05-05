#include <cassert>
#include <iostream>
#include <limits>

#include <cmath>


#include "../src/etops1.hpp"
#include "../src/etscalar.hpp"
#include "../src/etmatrix.hpp"
#include "../src/etops2.hpp"

#include "../src/array_template.hpp"
#include "../src/arrayops.hpp"
//#include "../src/tests_etarray.hpp"



using namespace std;
using std::cout;
using std::endl;


//const int THREADWORK = 2;

//const int THREADWORK = 1024;// << 5;
//const int THREADWORK = 1 << 12;

//const int maxThreads = 1 << 20;



// void checkResult(float *hostRef, float *gpuRef, const int N){
    
//     double epsilon = 1.0E-8;
//     bool match = 1;

//     for (int i = 0; i < N; i++){

//         if (abs(hostRef[i] - gpuRef[i]) > epsilon){
        
//             match = 0;
//             printf("Arrays do not match!\n");
//             printf("host %5.2f gpu %5.2f at current %d\n", hostRef[i],
//                    gpuRef[i], i);
//             break;
//         }//else{
//           //  printf("host %5.2f gpu %5.2f at current %d\n", hostRef[i],
//             //       gpuRef[i], i);
//         //}
//     }

//     if (match) printf("Arrays match.\n\n");

//     return;
// }


template <typename T>
void print (T const& c){

    int ncol = c.getncols();
    int nrow = c.getnrows();
    for (int i=0; i<nrow; ++i) {
        std::cout << '\n';
        for (int j=0; j<ncol; ++j) {
            std::cout << c(i,j) << ' ';
        }
    }
    printf("\ndone \n");
    //std::cout<< "\ndone \n" <<std::endl;
}

template <typename T, typename R>
void compare (T const& c, R const& z, 
                const std::string& input){
    // this string pass makes a copy ;-)

    int ncols = c.getncols();
    int nrows = c.getnrows();
    printf("checking simple  %s  \n",input.c_str());
    //std::cout << "checking simple " << input <<"\n";
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            assert( c(i,j) == z(i,j) );
        }
    }
    printf("simple %s OK \n\n",input.c_str());
    //std::cout << "simple " <<  input <<" OK \n\n";
}






//int etarraybasic (SolverManager& SM){
int etarraybasic (){



    int np = 10;

    //expression template matrices:
    Array<float> a(np,np), b(np,np), c(np,np), d(np,np);


    //compare with traditional overloaded matrices:
    Array2D<float> x(np,np), y(np,np), z(np,np), w(np,np);

    // initialize arrays with some values
    for (int i=0; i<np; ++i) {
         for (int j=0; j<np; ++j) {
            a(i,j) = static_cast <float>(i*(j+1));
            a(i,j) = 1.;//a(i,j) + 5.;
            b(i,j) = a(i,j)+a(i,j);
            c(i,j) = a(i,j)+b(i,j);

            x(i,j) = static_cast <float>(i*(j+1));
            x(i,j) = 1.;//x(i,j) + 5.;
            y(i,j) = x(i,j)+x(i,j);
            z(i,j) = x(i,j)+y(i,j);
         }
    }
    //...
    
    std::cout << "\n\n-----------------------------\n";
    std::cout << "\n\n-----------------------------\n";
    std::cout << "\n\n Welcome to the Test Folder! \n";
    std::cout << "\n\nTesting basic array class  \n";
    std::cout << "against  \n";
    std::cout << "Expression Template array class  \n";
    std::cout << "-----------------------------\n\n";
    std::cout << "\n\n-----------------------------\n";

    std::cout << "b: ";
    print(b);

    std::cout << "c: ";
    print(c);


    std::cout << "y: ";
    print(y);

    std::cout << "z: ";
    print(z);



    z = matmul(x,y);
    std::cout << "xfter z=matmul(x,y) ";
    print(z);
    d = 1.;

    c = matmul(a,b);
    std::cout << "after c=matmul(a,b) ";
    print(c);
    d = 1.;

    compare(c,z,"matmul");

    std::cout << "c = b + a: \n";
    std::cout.flush();
    c = b + a;
    std::cout << "c = b + a done \n";
    std::cout.flush();
    print(c);

    std::cout << "z = y + x: \n";
    std::cout.flush();
    z = y + x;
    std::cout << "z = y + x done \n";
    std::cout.flush();
    print(z);


    compare(c,z,"elementwise addition");

    
    a = 1.2f * a;  //expression template supports scalars
    std::cout << "1.2*a ";
    print(a);

    x = 1.2f * x;  //expression template supports scalars
    std::cout << "1.2*x ";
    print(x);


    compare(a,x,"scalar + matrix addition");

    std::cout << "a = c * a";
    a = c * a;
    print(a);
    std::cout << "a = c * a: ";
    print(a);


    std::cout << "x = z * x";
    x = z * x;
    print(x);
    std::cout << "x = z * x: ";
    print(x);

    compare(a,x,"elementwise matrix multiplication");


    a = 1.2f*a + a*b;
    std::cout << "1.2*a + a*b: ";
    print(a);


    x = 1.2f*x + x*y;
    std::cout << "1.2*x + x*y: ";
    print(x);

    compare(a,x,"scalar * matrix multiplication");

    a = b;
    std::cout << "after a = b: ";
    print(b);
    print(a);

    x = y;
    std::cout << "xfter x = y: ";
    print(y);
    print(x);

    compare(a,y,"crossed eualtiy");
    compare(b,x,"crossed equality");

    c = (a*b)+c;
    std::cout << "after c= (a*b)+c ";
    print(c);


    z = (x*y)+z;
    std::cout << "xfter z= (x*y)+z ";
    print(z);
    compare(c,z,"(mat * mat) + mat combined");


    c = (a+b)+c;
    std::cout << "after c= (a+b)+c ";
    print(c);
    
    z = (x+y)+z;
    std::cout << "xfter z= (x+y)+z ";
    print(z);

    compare(c,z,"(mat + mat) + mat combined");


    d=1.f;
    c = matmul(d,c);
    std::cout << "after c=matmul(d,c) ";
    print(c);

    w=1.f;
    z = matmul(w,z);
    std::cout << "xfter z=matmul(w,z) ";
    print(z);
    
    compare(c,z,"(matmul");
    



    std::cout << "b: ";
    print(b);
    //d = .5*b;
    compare(b,y,"b,y");

    //c = matmul(d,c);
    c = matmul((c+.5f*d),d); //combined & cached terms are supported
    std::cout << "after c=matmul((1.+b),c) ";
    print(c);

    
    std::cout << "y: ";
    print(y);
    //w = .5*y;
    compare(d,w,"d,w");

    //z = matmul(w,z);
    z = matmul((z+.5f*w),w); //compare with traditional overload
    std::cout << "xfter z=matmul((1.+y),z) ";
    print(z);


    compare(c,z,"compound mat + (scalar * mat )matmul mat");



    b=1;
    
    c=1.;
    c = matmul(b,c);
    std::cout << "after b=1.;c=1; c=matmul(b,c); c:";
    print(c);


    y=1;
    
    z=1.;
    z = matmul(y,z);
    std::cout << "xfter y=1.;z=1; z=matmul(y,z); z:";
    print(z);


    compare(c,z,"matmul");


    
    std::cout << "\n\n-----------------------------\n";
    std::cout << "\n\nTesting basic array class  \n";
    std::cout << "against  \n";
    std::cout << "Expression Template array class  \n";

    std::cout << "result:  \n";
    std::cout << "         all tests pass  \n";
    std::cout << "-----------------------------\n\n";


    
    return 0;

}