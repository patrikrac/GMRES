/*
Patrik Rác
Test File for the project in Grands systemes lineares.
The methods in this file test some part of the functionality of the code.
*/
#include "vector_op.hpp"
#include "dense_matrix.hpp"
#include "map_matrix.hpp"
#include "lu.hpp"

#include "solvers.hpp"

#include <cstdlib>
#include<iostream>
#include<fstream>

#define frand() ((double)rand() / RAND_MAX)*2.0 - 1.0

bool testLU();

bool testGMRES();


int main(int argc, char *argv[])
{
    /*Call the different test funcitons and check return value.*/
    std::cout << "Testing..." << std::endl;
    std::cout << "Running test for the complex LU solver..." << std::endl;
    std::cout << (testLU() ? "Success." : "Failure.") << std::endl;

    std::cout << "Running test for the complex GMRES solver..." << std::endl;
    std::cout << (testGMRES() ? "Success." : "Failure.") << std::endl;
}

bool testLU()
{
    int size = 512;
    srand((unsigned) time(NULL)); 
    DenseMatrix<Cplx> _A_TEST1(size,size);
    std::vector<Cplx> _b_TEST1(size);
    std::vector<Cplx> _x_TEST1(size);
    for(int j = 0; j < size; j++)
    {
        for(int i = 0; i < size; i++)
        {
            _A_TEST1(i,j) = Cplx(frand(), frand());
        }
        _x_TEST1[j] = Cplx(frand(), frand());
    }
    _b_TEST1 = _A_TEST1*_x_TEST1;

    LUSolver<Cplx> _lu_TEST1(_A_TEST1);
    std::vector<Cplx > _res_TEST1 = _lu_TEST1.Solve(_b_TEST1);
    for(int i = 0; i < size; i++)
    {
        //Check the difference of the exact solution with the copmuted one with a fixed threshhold
        if(std::abs(std::abs(_res_TEST1[i]) - std::abs(_x_TEST1[i])) > 1e-3) return false;
    }

    return true;
}

bool testGMRES()
{
    /*Test the GMRES Solver on a random small dense problem wiht restart = size*/
    int size = 64;
    srand((unsigned) time(NULL)); 
    MapMatrix<Cplx> _A_TEST1(size,size);
    std::vector<Cplx> _b_TEST1(size);
    std::vector<Cplx> _x_TEST1(size);
    for(int j = 0; j < size; j++)
    {
        for(int i = 0; i < size; i++)
        {
            _A_TEST1.insert(i,j,Cplx(frand(), frand())); 
        }
        _x_TEST1[j] = Cplx(frand(), frand());
    }
    _b_TEST1 = _A_TEST1*_x_TEST1;

    //Copmute solution using GMRES
    std::vector<Cplx > _res_TEST1(_A_TEST1.NbCol(), 1.0);
    GMResSolve(_A_TEST1, _res_TEST1, _b_TEST1, size, 1e-4, 1000);
    for(int i = 0; i < size; i++)
    {
        //Check the difference of the exact solution with the copmuted one with a fixed threshhold
        if(std::abs(std::abs(_res_TEST1[i]) - std::abs(_x_TEST1[i])) > 1e-3) return false;
    }

    /*Test the GMRES on the matrix 1 of the matrix archive with a random solution*/
    MapMatrix<Cplx> _A_TEST2 = LoadMapMatrix<Cplx>("../matrix-archive/matrix_1.txt");
    std::vector<Cplx> _b_TEST2(_A_TEST2.NbCol());
    std::vector<Cplx> _x_TEST2(_A_TEST2.NbRow());
    for(int i = 0; i < _x_TEST2.size(); i++)
    {
        _x_TEST2[i] = Cplx(frand(), frand());
    }

    _b_TEST2 = _A_TEST2*_x_TEST2;

     //Copmute solution using GMRES
    std::vector<Cplx > _res_TEST2(_A_TEST2.NbCol(), 1.0);
    GMResSolve(_A_TEST2, _res_TEST2, _b_TEST2, 80, 1e-9, 1e4);
    for(int i = 0; i < size; i++)
    {
        //Check the difference of the exact solution with the copmuted one with a fixed threshhold
        if(std::abs(std::abs(_res_TEST2[i]) - std::abs(_x_TEST2[i])) > 1e-3) return false;
    }

    return true;
}
