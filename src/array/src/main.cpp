
#include <cassert>
#include <iostream>
#include <limits>

#include "array_template.hpp"
#include "SolverManager.hpp"
#include "input.hpp"
#include "geometry.hpp"

#include "tests_array.hpp"
#include "tests_etarray.hpp"


using namespace std;





int main (int argc, char **argv) {




    SolverManager SM = SolverManager();
    SM.driver(argc, argv);


    char type;
    do
    {
        cout << "Run some (quick!) tests? (y/n)" << endl;
        cin >> type;
    }
    while( !cin.fail() && type!='y' && type!='n' );

    if ( type =='y' ) {
        //array basic tests
        arraybasic();
        //expression tempalte (et) array basic tests
        etarraybasic();
        //etarraybasic(SM);
    }


    return 0;
}