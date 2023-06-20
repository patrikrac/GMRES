/*
Patrik Rác
Grands systemes lineares Projet
Main file printing the convergence results of the given matrices to a file compatible with gnuplot
*/
#include "map_matrix.hpp"
#include "solvers.hpp"

#include<iostream>
#include<chrono>
#include<complex>

typedef std::complex<double> Cplx;

/*Funciton benchmarking the convergence of the GMRES solver, takes the restart parameter as an argument*/
void convergenceRun(int restart)
{
    /*Setting parameters for the data input and output files.*/
    std::string in_A("../matrix-archive/matrix_");
    std::string in_b("../matrix-archive/rhs_");
    
    std::string out_conv("convergence_");

    for (int i = 1; i <= 5; i++)
    {
        std::string mat_A = in_A + std::to_string(i) + std::string(".txt");
        std::string vec_b = in_b + std::to_string(i) + std::string(".txt");
        std::string fn = out_conv + std::to_string(i) +std::string("_") + std::to_string(restart) + std::string(".txt");

        std::ofstream out_file(fn);

        MapMatrix<Cplx> A = LoadMapMatrix<Cplx>(mat_A);
        std::vector<Cplx> b = LoadVector<Cplx>(vec_b);

        std::vector<Cplx> x(A.NbCol(), 1.0);

        std::cout << "Computing equation " << i << std::endl;

        auto begin = std::chrono::high_resolution_clock::now();
        GMResSolve(A, x, b, restart, 1e-6, 1e4, out_file);
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::cout << "Done in " << elapsed.count()*1e-9 << "s" << std::endl;

        out_file.close();
    }
}

int main(int argc, char *argv[])
{
    /*Some runs of the convergence benchmarks with different values of restart*/
    convergenceRun(20);
    convergenceRun(40);
    convergenceRun(80);
    return 0;
}
