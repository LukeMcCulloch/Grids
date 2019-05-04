//
//=================================
// include guard
#ifndef __SOLVERMANAGER_INCLUDED__
#define __SOLVERMANAGER_INCLUDED__

//#include "input.hpp"
//#include "geometry.hpp"


#include "tests_etarray.hpp"


using namespace std;







class SolverManager{

    public:

        void input_discrete(int argc, char **argv);
        void driver(int argc, char **argv);
        void plate();
        void startMessage();
        

    private:

        int NI; // Number of x pts
        int NJ; // Number of y pts
        int nt; // number of time steps

        float  lx;  // Length of Plate in X
        float  ly;  // Length of Plate in Y
        float  Xo;  // 00 x
        float  Yo;  // 00 y
        float  C;   // courant number
        float  Re;  // reynolds num
};


void SolverManager::driver(int argc, char *argv[]){

    startMessage();


    input_discrete(argc, argv);

    cout << "Hello plate?" << endl;
    plate();

    
}

void SolverManager::startMessage(){
    cout << "-------------------------------------------\n";
    cout << "Startup:" << endl;
    cout << "matrix algebra test driver" << endl;
    cout << "dummy solver program, and then test matrix operations" << endl;
    cout << "" << endl;
    cout << "-------------------------------------------\n" ;
}



#endif //__TESTS_ARRAY_INCLUDED__