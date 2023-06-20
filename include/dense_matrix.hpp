/*
Patrik Rác
Grands systemes lineares
Implementation of the dense matrix class.
*/
#pragma once
#include"map_matrix.hpp"

#include<cassert>
#include<iostream>
#include<vector>

/*Forward declarations of friend function of templated dense matrix*/
template<typename value_t>
class DenseMatrix;

template<typename value_t>
std::ostream &operator<<(std::ostream &, DenseMatrix<value_t> &);

template<typename value_t>
DenseMatrix<value_t> operator+(const DenseMatrix<value_t> &, const DenseMatrix<value_t> &);

template<typename value_t>
DenseMatrix<value_t> operator*(const DenseMatrix<value_t> &, const DenseMatrix<value_t> &);

template<typename value_t>
DenseMatrix<value_t> operator*(const value_t &, const DenseMatrix<value_t> &);

template<typename value_t>
std::vector<value_t> operator*(const DenseMatrix<value_t> &, const std::vector<value_t> &);

/*Class DenseMatrix modelling a densly populated matrix with elements of type value_t*/
template<typename value_t>
class DenseMatrix
{
    int nr;
    int nc;
    std::vector<value_t> data;

public:
    DenseMatrix(const int &nr, const int &nc) : nr(nr), nc(nc), data(nr*nc)
    {}

    DenseMatrix(const DenseMatrix &);

    DenseMatrix &operator=(const DenseMatrix &);
    value_t &operator()(const int, const int);
    const value_t &operator()(const int, const int) const;

    friend std::ostream &operator<<<value_t>(std::ostream &, DenseMatrix &);

    friend DenseMatrix operator+<value_t>(const DenseMatrix &, const DenseMatrix &);

    DenseMatrix &operator+=(const DenseMatrix &);

    friend DenseMatrix operator*<value_t>(const DenseMatrix &, const DenseMatrix &);

    friend DenseMatrix operator*<value_t>(const value_t &, const DenseMatrix &);

    friend std::vector<value_t> operator*<value_t>(const DenseMatrix &, const std::vector<value_t> &);

    DenseMatrix &operator*=(const DenseMatrix &);

    DenseMatrix &operator*=(const value_t &);

    bool operator==(const DenseMatrix &);

    int NbRow() {return nr;}

    int NbCol() {return nc;}

    friend int NbRow(const DenseMatrix &m){return m.nr;}

    friend int NbCol(const DenseMatrix &m) {return m.nc;}
};

template<typename value_t>
DenseMatrix<value_t>::DenseMatrix(const DenseMatrix<value_t> &m)
{
    nc = m.nc;
    nr = m.nr;
    data = m.data;
}

template<typename value_t>
DenseMatrix<value_t> &DenseMatrix<value_t>::operator=(const DenseMatrix<value_t> &m)
{
    nc = m.nc;
    nr = m.nr;
    data = m.data;

    return *this;
}

template<typename value_t>
value_t &DenseMatrix<value_t>::operator()(const int i, const int j)
{
    assert(i < nr && j < nc);
    return data[i*nc + j];
}

template<typename value_t>
const value_t &DenseMatrix<value_t>::operator()(const int i, const int j) const 
{
    assert(i < nr && j < nc);
    return data[i*nc + j];
}

template<typename value_t>
std::ostream &operator<<(std::ostream &stream, DenseMatrix<value_t> &m)
{
    for(int i = 0; i < m.nr; i++)
    {
        for(int j = 0; j < m.nc; j++)
        {
            stream << m(i,j) << "\t";
        }
        stream << std::endl;
    }
    return stream;
}

template<typename value_t>
DenseMatrix<value_t> operator+(const DenseMatrix<value_t> &m1, const DenseMatrix<value_t> &m2)
{
    assert(m1.nc==m1.nc && m2.nr==m2.nr);
    DenseMatrix res(m1.nr, m1.nc);

    for (int i = 0; i < m1.nc*m1.nr; i++)
    {
        res.data[i] = m1.data[i]+m2.data[i];
    }

    return res;
}

template<typename value_t>
DenseMatrix<value_t> &DenseMatrix<value_t>::operator+=(const DenseMatrix<value_t> &m)
{
    assert(nc==m.nc && nr==m.nr);

    for (int i = 0; i < nc*nr; i++)
    {
        data[i] += m.data[i];
    }
    return *this;
}

template<typename value_t>
DenseMatrix<value_t> operator*(const DenseMatrix<value_t> &m1, const DenseMatrix<value_t> &m2)
{
    assert(m1.nc==m2.nr);
    DenseMatrix<value_t> res(m1.nr, m2.nc);

    for(int i = 0; i < res.nr; i++)
    {
        for(int j = 0; j < res.nc; j++)
        {
            for(int k = 0; k < m1.nc; k++)
            {
                res(i,j) += m1(i,k)*m2(k,j);
            }
        }
    }

    return res;
}

template<typename value_t>
DenseMatrix<value_t> operator*(const value_t &a, const DenseMatrix<value_t> &m)
{
    DenseMatrix res(m);

    for (int i=0; i < res.nc*res.nr; i++)
    {
        res.data[i] *= a;
    }

    return res;
}

template<typename value_t>
std::vector<value_t> operator*(const DenseMatrix<value_t> &m, const std::vector<value_t> &v)
{
    assert(v.size() == m.nc);

    std::vector<value_t> res(m.nr);

    for(int i = 0; i < m.nr; i++)
    {
        for(int k = 0; k < m.nc; k++)
        {
            res[i] += m(i,k)*v[k]; 
        }
    }

    return res;
}

template<typename value_t>
DenseMatrix<value_t> &DenseMatrix<value_t>::operator*=(const DenseMatrix<value_t> &m)
{
    assert(nc==m.nr);

    std::vector<value_t> tmp(nr*m.nc);

    for(int i = 0; i < nr; i++)
    {
        for(int j = 0; j < m.nc; j++)
        {
            for(int k = 0; k < nc; k++)
            {
                tmp[i*m.nc + j] += this->operator()(i,k)*m(k,j);
            }
        }
    }

    //Adjust the own members
    nc = m.nc;
    data = tmp;

    return *this;
}

template<typename value_t>
DenseMatrix<value_t> &DenseMatrix<value_t>::operator*=(const value_t &a)
{
    for (int i=0; i < nc*nr; i++)
    {
        data[i] *= a;
    }

    return *this;
}

template<typename value_t>
bool DenseMatrix<value_t>::operator==(const DenseMatrix &m)
{
    assert(this->nr == m.nr && this->nc == m.nc);
    for(int i = 0; i < nr*nc; i++)
    {
        if(data[i] != m.data[i]) return false;
    }
    return true;
}