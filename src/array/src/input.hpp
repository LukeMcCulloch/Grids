//=================================
// include guard
#ifndef __INPUT_INCLUDED__
#define __INPUT_INCLUDED__

// basic file operations
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
//#include <iomanip>

//#include "input.hpp"
#include "SolverManager.hpp"

using namespace std;


void SolverManager::input_discrete(int argc, char *argv[]){ 
    // user input:
    cout << "input_discrete " <<endl;
    ifstream infile;
    //streampos begin,end;
    string infilename;

    cout << "Please enter the input file name> " << flush;
    cout << "Default is fifi.dat" << endl;
    cout << "(right now this does nothing much.  just press enter)" << endl;
    while (true){

        getline( cin, infilename );
        if (infilename.length()>1){

            cout << "Opening " << infilename <<endl;
        }
        else{

            infilename = "fifi.dat";
            cout << "Opening " << infilename <<endl;
        }
        infile.open( infilename.c_str() );
        if (infile) break;
        cout << "Invalid file. Please enter a valid input file name> " << flush;
    }



    string line;

    if (infile.is_open()){
        // begin = infile.tellg();
        // infile.seekg (0, ios::end);
        // end = infile.tellg();
        // infile.seekg (0, ios::beg);
        getline (infile,line);
            cout << "title: " << line << '\n';
        //getline (infile,line);
            //cout << line << '\n';
            infile >> this->Xo >> this->Yo;
            infile >> this->lx >> this->ly;
            infile >> this->NI >> this->NJ >> this->nt;
            infile >> this->C;
            infile >> this->Re;
            // cout << "Xo = " << this->Xo << " Yo = "<< this->Yo << '\n';
            // cout << "lx = " << this->lx << " ly = "<< this->ly << '\n';
            // cout << "NI = " << this->NI << " NJ = "<< this->NJ << " nt = "<< this->nt<< '\n';
            // cout << "C = " << this->C << '\n';
            // cout << "Re = " << this->Re << '\n';
        while ( getline (infile,line) ){

            cout << line << '\n';
        }
        
        infile.close();
    }

    else{ 
      cout << "Unable to open file"; 
    }
}


#endif //__MYCLASS_H_INCLUDED__


/*
https://stackoverflow.com/questions/25274312/
        defining-c-functions-inside-header-files

If you want to use a function in multiple source files 
(or rather, translation units), then you place 
a function declaration (i.e. a function prototype) 
in the header file, and the definition in one source file.


*/