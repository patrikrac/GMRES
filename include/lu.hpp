/*
Patrik Rác
Grands systemes lineares projet
Implementation of a templated LU solver class.
The class is designed to work with copmlex valued types.
The main class uses partial pivoting to improve accuracy.
*/
#pragma once
#include "dense_matrix.hpp"
#include <vector>
#include <cassert>
#include <iostream>
#include <complex>


/*Function that turns a Dense matrix into its LU factorisation without pivoting*/
template<typename value_t>
void LUFactorize(DenseMatrix<value_t> &);

/*
Function that turns a Dense matrix into its LU factorisation with partial pivoting
Permutation is stored in  a vector
*/
template<typename value_t>
void LUFactorize(DenseMatrix<value_t> &, std::vector<int> &);

/*Forward declaration of the LU solver class*/
template<typename value_t>
class LUSolver;

template<typename value_t>
std::ostream &operator<<(std::ostream &o, const LUSolver<value_t> &lu);

template<typename value_t>
class LUSolver
{
    int nr, nc;
    std::vector<int> pivot;
    DenseMatrix<value_t> LU;

public:
    LUSolver(const DenseMatrix<value_t> &);
    LUSolver(const LUSolver &lu):nr(lu.nr), nc(lu.nc), pivot(lu.pivot), LU(lu.LU)
    {}

    LUSolver &operator=(const LUSolver &lu)
    {
        nr = lu.nr;
        nc = lu.nc;
        pivot = lu.pivot;
        LU = lu.LU;

        return *this;
    }

    const value_t &operator()(const int &j, const int &k) const;

    friend std::ostream &operator<<<value_t>(std::ostream &o, const LUSolver &lu);

    std::vector<value_t> Solve(const std::vector<value_t>& b);

    friend int NbRow(const LUSolver &lu)  {return lu.nr;}
    friend int NbCol(LUSolver &lu) {return lu.nc;}
};

template<typename value_t>
void LUFactorize(DenseMatrix<value_t> &mat)
{
    int n = NbCol(mat);
    value_t pivot;

    for(int i = 0; i < n; i++)
    {
        pivot = mat(i,i);
        for (int j = i+1; j < n; j++)
        {
            /*Update the lower part of the matrix*/
            mat(j,i) /= pivot; 

            /*Update the upper part of the matrix iteratively*/
            for(int k = i+1; k < n; k++)
            {
                mat(j,k) -=  mat(j,i)*mat(i,k);
            }
        }
    }
}

/*Helper function that swaps the rows i and k of a matrix mat*/
template<typename value_t>
static void swapRows(DenseMatrix<value_t> &mat, const int &j, const int &k)
{
    value_t tmp;
    for (int i = 0; i < NbCol(mat); i++)
    {
        tmp = mat(j,i);
        mat(j,i) = mat(k,i);
        mat(k,i) = tmp;
    }
}

template<typename value_t>
void LUFactorize(DenseMatrix<value_t> &mat, std::vector<int> &pivot)
{
    int n = NbCol(mat);

    for(int i = 0; i < n; i++)
    {
        /*Find pivot element for current column*/
        value_t p = mat(i,i);
        int pivotIdx = i;
        for(int k = i+1; k < NbRow(mat); k++)
        {
            if(std::abs(p) < std::abs(mat(k,i))) 
            {
                p = mat(k,i);
                pivotIdx = k;
            }
        }

        /*Swap the  rows if pivoting is required*/
        if (i != pivotIdx)
        {
            int tmp = pivot[i];
            pivot[i] = pivot[pivotIdx];
            pivot[pivotIdx] = tmp;

            swapRows(mat, i, pivotIdx);
        }
        
        /*Execute regular LU Factorisation*/
        for (int j = i+1; j < n; j++)
        {
            /*Update the lower part of the matrix*/
            mat(j,i) = mat(j,i)/p; 

            /*Update the upper part of the matrix iteratively*/
            for(int k = i+1; k < n; k++)
            {
                mat(j,k) -= mat(j,i)*mat(i,k);
            }
        }
    }
}

/*Implementation of the LUSolver class functions*/
template<typename value_t>
LUSolver<value_t>::LUSolver(const DenseMatrix<value_t> &mat) : nr(NbRow(mat)), nc(NbCol(mat)), pivot(nr), LU(mat)
{
    for(int i = 0; i < nr; i++) pivot[i] = i;
    LUFactorize(LU, pivot);
}

template<typename value_t>
const value_t &LUSolver<value_t>::operator()(const int &j, const int &k) const
{
    assert(j<nr && k<nc);
    return LU(j,k);
}

template<typename value_t>
std::ostream &operator<<(std::ostream &o, const LUSolver<value_t> &lu)
{
    o << "L =";
    for (int i = 0; i < lu.nr; i++)
    {
        o << "\t";
        for(int j = 0; j < i; j++)
        {
            o << lu(i, j) << "\t";
        }   
        o << 1 << "\t";
        for(int j = i+1; j < lu.nr; j++)
        {
            o << 0 << "\t";
        }  
        o << "\n";
    }
    o << std::endl;

    o << "U =";
     for (int i = 0; i < lu.nr; i++)
    {
        o << "\t";
        for(int j = 0; j < i; j++)
        {
            o << 0 << "\t";
        }   
        for(int j = i; j < lu.nr; j++)
        {
            o << lu(i, j) << "\t";
        }   
        o << "\n";
    }
    o << std::endl;
    o << "P =\t";
    for(int i = 0; i < lu.nr; i++)
    {
        o << lu.pivot[i] << "\t";
    }
    o << std::endl;

    return o;
}

/*Solve the system PLUx = b using the created partial pivoted LU factorisation*/
template<typename value_t>
std::vector<value_t> LUSolver<value_t>::Solve(const std::vector<value_t> &b)
{
    assert(b.size() == nc);


    std::vector<value_t> sol(b.size());

    /*Apply P to the right hand side vector b*/
    for(int i = 0; i < b.size(); i++) sol[i] = b[pivot[i]];

    /*Solve the System Ly = b*/
    for(int i = 0; i < nc; i++)
    {
        for (int j = i+1; j < nc; j++)
        {
            sol[j] -= sol[i]*this->LU(j,i); 
        } 
    }

    /*Solve the system Ux = y*/
    for(int i = nc-1; i >= 0; i--)
    {
        sol[i] /= this->LU(i,i);
        for(int j = i-1; j >= 0; j--)
        {
            sol[j] -= sol[i]*this->LU(j,i);
        }
    }

    return sol;
}