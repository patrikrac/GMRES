/*
Patrik Rác
Grands systemes lineares projet
Implentation of the solver functions. Especially the NormalSolve and the GMRES-solver.
*/

#include "solvers.hpp"
#include "vector_space.hpp"

/*Computes a hermitian transpose of a given complex dense matrix A*/
static DenseMatrix<Cplx> hermitianTranspose(const DenseMatrix<Cplx> &A)
{
    DenseMatrix<Cplx> Aht(NbCol(A), NbRow(A));
    for(int i = 0; i < NbRow(A); i++)
    {
        for(int j = 0; j < NbCol(A); j++)
        {
            Aht(j,i) = std::conj(A(i,j));
        }
    }
    return Aht;
}

/*Solves the least squares problem ||Ax - b||*/
void NormalSolve(const DenseMatrix<Cplx> &A, std::vector<Cplx> &x, const std::vector<Cplx> &b)
{
    /*Compute the hermitian transpose of A*/
    DenseMatrix<Cplx> Aht = hermitianTranspose(A);
    /*Compute the new right hand side A*b*/
    std::vector<Cplx> c = Aht*b; 
    /*Solve the linear system A*Ax = A*b using the LU solver*/
    LUSolver<Cplx> lu(Aht*A);
    x = lu.Solve(c);    
}

/*Copmute Arnoldi computing V and H following a given vector v_1 in the vector space V that is assumed to be normalized*/
static void Arnoldi(const MapMatrix<Cplx> &A, VectorSpace<Cplx> &V, DenseMatrix<Cplx> &H, const int &m)
{
    std::vector<Cplx> Av;
    std::vector<Cplx> omega;
     for (int j = 0; j < m; j++)
    {
        Av = A*V[j];
        omega = Av;
        for(int i = 0; i <= j; i++) 
        {
            H(i,j) = (Av,V[i]);
            omega -= H(i,j) * V[i];
        }
        
        H(j+1,j) = Norm(omega);

        if(H(j+1,j) == 0.0) break;
        
        V[j+1] = (1.0/H(j+1,j))*omega;
    }
}

/*Restarted GMRES funciton that solves Ax=b iteratively */
void GMResSolve(const MapMatrix<Cplx> &A, std::vector<Cplx> &x, const std::vector<Cplx> &b, int restart, double tol, int maxit)
{
    /*Initialization*/
    VectorSpace<Cplx> V(restart+1, x.size());
    std::vector<Cplx> y(restart);
    std::vector<Cplx> e1(restart+1, 0.0);
    e1[0] = 1.0;

    DenseMatrix<Cplx> H(restart+1, restart);
    
    for(int iter = 0; iter < maxit; iter++)
    {
        /*Compute the residual and its norm*/
        std::vector<Cplx> r = b - A*x;
        double beta = Norm(r);
        V[0] = (1.0/beta)*r;
        
        /*Check stopping condition*/
        if(beta/Norm(b) < tol) break;

        /*Generate the Arnoldi basis and the matrix H using Arnoldi*/
        Arnoldi(A, V, H, restart);

        /*Compute the solution to the least squares problem ||beta*e_1 - H*y|| => y_m and x_m = x_0 + V_m*y_m*/
        NormalSolve(H, y, -beta*e1);
        
        /*Update the solution vector by x_m = x_0 + V_m*y_m*/
        x += V*y;
    }
}

/*Overload of the GMRES solves that plots convergence history to a file specified by log*/
void GMResSolve(const MapMatrix<Cplx> &A, std::vector<Cplx> &x, const std::vector<Cplx> &b, int restart, double tol, int maxit, std::ofstream &log)
{
    /*Initialization*/
    VectorSpace<Cplx> V(restart+1, x.size());
    std::vector<Cplx> y(restart);
    std::vector<Cplx> e1(restart+1, 0.0);
    e1[0] = 1.0;

    double res;

    DenseMatrix<Cplx> H(restart+1, restart);
    
    for(int iter = 0; iter < maxit; iter++)
    {
        /*Compute the residual and its norm*/
        std::vector<Cplx> r = b - A*x;
        double beta = Norm(r);
        V[0] = (1.0/beta)*r;
        
        /*Log the quadratic norm of the relative residual*/
        log << iter << "\t" << (res = beta/Norm(b)) << std::endl;

         /*Check stopping condition*/
        if(res < tol) break;
        
        /*Generate the Arnoldi basis and the matrix H using Arnoldi*/
        Arnoldi(A, V, H, restart);

        /*Compute the solution to the least squares problem ||beta*e_1 - H*y|| => y_m and x_m = x_0 + V_m*y_m*/
        NormalSolve(H, y, -beta*e1);

        /*Update the solution vector by x_m = x_0 + V_m*y_m*/
        x += V*y;
    }

}