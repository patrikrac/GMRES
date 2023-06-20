/*
Patrik Rác
Grands systemes lineares
Implementaion of Question 2 and Question 3 of the Project sheet
*/

#pragma once
#include "vector_op.hpp"
#include "dense_matrix.hpp"
#include "map_matrix.hpp"
#include "lu.hpp"

#include <vector>
#include<complex>

typedef std::complex<double> Cplx;

/*Applies method of the normal equation to solve the least squares problem ||Ax - b||*/
void NormalSolve(const DenseMatrix<Cplx> &A, std::vector<Cplx> &x, const std::vector<Cplx> &b);

/*Models a restarted GMRes method for solving the linear system Ax=b*/
void GMResSolve(const MapMatrix<Cplx> &A, std::vector<Cplx> &x, const std::vector<Cplx> &b, int restart, double tol, int maxit);

/*Overload of the GMResSolve with a log output-stream to be handeled by the user*/
void GMResSolve(const MapMatrix<Cplx> &A, std::vector<Cplx> &x, const std::vector<Cplx> &b, int restart, double tol, int maxit, std::ofstream &log);