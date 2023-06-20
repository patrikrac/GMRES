/*
Patrik Rác
Grands systemes lineares projet
Implementation of templated vector operations.
*/
#pragma once
#include<vector>
#include<math.h>
#include<cassert>
#include<complex>

#include<iostream>

#include<fstream>

typedef std::complex<double> Cplx;

/*Forward declaration of templated vector operations*/

template <typename value_t>
std::vector<value_t> &operator+=(std::vector<value_t> &, const std::vector<value_t> &);

template <typename value_t>
std::vector<value_t> &operator-=(std::vector<value_t> &, const std::vector<value_t> &);

template <typename value_t>
std::vector<value_t> operator+(const std::vector<value_t> &, const std::vector<value_t> &);

template <typename value_t>
std::vector<value_t> operator-(const std::vector<value_t> &, const std::vector<value_t> &);

template <typename value_t>
value_t operator,(const std::vector<value_t> &, const std::vector<value_t> &);

template <>
inline Cplx operator,(const std::vector<Cplx> &, const std::vector<Cplx> &);

template<typename value_t>
double Norm(const std::vector<value_t> &);

template<>
inline double Norm(const std::vector<Cplx> &);

template <typename value_t>
std::vector<value_t> &operator*=(std::vector<value_t> &, const value_t &);

template <typename mult_value_t, typename value_t>
std::vector<value_t> operator*(const mult_value_t &, const std::vector<value_t> &);

template <typename value_t>
std::vector<value_t> LoadVector(const std::string& filename);

template<>
inline std::vector<Cplx> LoadVector(const std::string& filename);

template<typename value_t>
void Print(std::vector<value_t> &v);


/*End of forward declaration*/


template <typename value_t>
std::vector<value_t> &operator+=(std::vector<value_t> &v1, const std::vector<value_t> &v2)
{
    assert(v1.size()==v2.size());

    for (int i = 0; i < v1.size(); i++)
    {
        v1[i] += v2[i];
    }

    return v1;
}

template <typename value_t>
std::vector<value_t> &operator-=(std::vector<value_t> &v1, const std::vector<value_t> &v2)
{
    assert(v1.size()==v2.size());

    for (int i = 0; i < v1.size(); i++)
    {
        v1[i] -= v2[i];
    }

    return v1;
}

template <typename value_t>
std::vector<value_t> operator+(const std::vector<value_t> &v1, const std::vector<value_t> &v2)
{
    assert(v1.size()==v2.size());
    std::vector<value_t> res(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        res[i] = v1[i] + v2[i];
    }

    return res;
}

template <typename value_t>
std::vector<value_t> operator-(const std::vector<value_t> &v1, const std::vector<value_t> &v2)
{
    assert(v1.size()==v2.size());
    std::vector<value_t> res(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        res[i] = v1[i] + v2[i];
    }

    return res;
}

template <typename value_t>
value_t operator,(const std::vector<value_t> &v1, const std::vector<value_t> &v2)
{
    assert(v1.size()==v2.size());
    value_t sum = (value_t) 0.0;

    for (int i = 0; i < v1.size(); i++)
    {
        sum += v1[i]*v2[i];
    }

    return sum;
}

template <>
inline Cplx operator,(const std::vector<Cplx> &v1, const std::vector<Cplx> &v2)
{
    assert(v1.size()==v2.size());
    Cplx sum =  0.0;

    for (int i = 0; i < v1.size(); i++)
    {
        sum += v1[i]*std::conj(v2[i]);
    }

    return sum;
}


template<typename value_t>
double Norm(const std::vector<value_t> &v)
{
    double sum = 0.0;
    for(int i = 0; i < v.size(); i++)
    {
        sum += v[i]*v[i];
    }

    return sqrt(sum);
}

template<>
inline double Norm(const std::vector<Cplx> &v)
{
    Cplx sum = 0.0;
    for(int i = 0; i < v.size(); i++)
    {
        sum += v[i]*std::conj(v[i]);
    }

    return sqrt(std::real(sum));
}


template <typename value_t>
std::vector<value_t> &operator*=(std::vector<value_t> &v, const double &a)
{
    for(int i = 0; i < v.size(); i++)
    {
        v[i] *= a;
    }

    return v;
}

template <typename mult_value_t, typename value_t>
std::vector<value_t> operator*(const mult_value_t &a, const std::vector<value_t> &v)
{
    std::vector<value_t> res(v.size());
    for(int i = 0; i < v.size(); i++)
    {
        res[i] = a*v[i];
    }

    return res;
}

template <typename value_t>
std::vector<value_t> LoadVector(const std::string& filename)
{
    std::ifstream file(filename);

    std::vector<value_t> v;

    value_t elem;

    while(file >> elem)
    {
        v.push_back(elem);
    }

    file.close();
    return v;
}


template <>
inline std::vector<Cplx> LoadVector(const std::string& filename)
{
    std::ifstream file(filename);

    std::vector<Cplx> v;

    double elem_real, elem_imag;

    while(file >> elem_real >> elem_imag)
    {
        v.push_back(Cplx(elem_real, elem_imag));
    }

    file.close();
    return v;
}

template<typename value_t>
void Print(std::vector<value_t> &v)
{
    for(int i = 0; i < v.size(); i++)
    {
        std::cout << v[i] << "\t";
    }
    std::cout << std::endl;
}
