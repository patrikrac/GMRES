/*
Patrik Rác
Grands systemes lineares projet
Simple additional class that manages a vector space
*/
#pragma once
#include<vector>
#include<cassert>

template<typename value_t>
class VectorSpace;

template<typename value_t>
class VectorSpace
{
    /*The class manages n vectors of dimension d*/
    int n, d;

    std::vector<std::vector<value_t>> space;

public:
    VectorSpace(const int n, const int d): n(n), d(d), space(n)
    {}

    std::vector<value_t> &operator[](const int &i);
    const std::vector<value_t> &operator[](const int &i) const;

    std::vector<value_t> operator*(std::vector<value_t> &r);

    int size() {return n;}
};

template<typename value_t>
std::vector<value_t> &VectorSpace<value_t>::operator[](const int &i)
{
    assert((0 <= i) && (i < n));
    return space[i];
}

template<typename value_t>
const std::vector<value_t> &VectorSpace<value_t>::operator[](const int &i) const
{
    assert((0 <= i) && (i < n));
    return space[i];
}

/*Multiply vector space with a given vector r*/
template<typename value_t>
std::vector<value_t> VectorSpace<value_t>::operator*(std::vector<value_t> &r)
{
    assert(r.size() <= n);

    std::vector<value_t> res(d, 0.0);

    for(int j = 0; j < r.size(); j++)
    {
        for(int i = 0; i < d; i++)
        {
            res[i] += space[j][i] * r[j];
        }
    }
    return res;
}
