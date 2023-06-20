/*
Patrik Rác
Grands systems lineares projet
Implementation of the templated MapMatrix class.
*/
#pragma once

#include<tuple>
#include<map>
#include<vector>
#include<iostream>
#include<cassert>
#include<complex>

#include<string>
#include<fstream>

typedef std::complex<double> Cplx;

/*Forward declarations of the MapMatrix class*/
template <typename value_t>
class MapMatrix;

/*Declaring templated friend functions*/
template <typename value_t>
MapMatrix<value_t> operator*(const MapMatrix<value_t> &, const MapMatrix<value_t> &);

template <typename value_t>
MapMatrix<value_t> operator*(const value_t &, const MapMatrix<value_t> &);

template <typename value_t>
std::vector<value_t> operator*(const MapMatrix<value_t> &, const std::vector<value_t> &);

template <typename value_t> 
std::ostream &operator<<(std::ostream &, MapMatrix<value_t> &);

template <typename value_t>
MapMatrix<value_t> operator+(const MapMatrix<value_t> &, const MapMatrix<value_t> &);

template <typename value_t>
MapMatrix<value_t> LoadMapMatrix(const std::string& filename);

template <>
inline MapMatrix<Cplx> LoadMapMatrix(const std::string& filename);

template <typename value_t>
void Write(const std::string &filename, MapMatrix<value_t>& m);

template <>
inline void Write(const std::string &filename, MapMatrix<Cplx>& m);


/*Declarations of the Map Matrix class and its member functions*/
template<typename value_t>
class MapMatrix
{
    int nr;
    int nc;
    typedef std::tuple<int,int> NxN;
    std::map<NxN, value_t> data;

public:
    MapMatrix(const int &nr, const int &nc) : nr(nr), nc(nc) 
    {}

    MapMatrix(const MapMatrix &m): nr(m.nr), nc(m.nc), data(m.data) 
    {}

    MapMatrix &operator=(const MapMatrix &);

    void insert(const int &, const int &, const value_t &);

    friend std::ostream &operator<<<value_t>(std::ostream &, MapMatrix &);

    friend MapMatrix operator+<value_t>(const MapMatrix &, const MapMatrix &);

    MapMatrix &operator+=(const MapMatrix &);

    friend MapMatrix operator*<value_t>(const MapMatrix &, const MapMatrix &);

    friend MapMatrix operator*<value_t>(const value_t &, const MapMatrix &);

    friend std::vector<value_t> operator*<value_t>(const MapMatrix &, const std::vector<value_t> &);

    MapMatrix &operator*=(const MapMatrix &);

    MapMatrix &operator*=(const double &);

    int NbRow() const {return nr;} 

    int NbCol() const {return nc;}

    friend MapMatrix LoadMapMatrix<value_t>(const std::string& filename);
    friend void Write<value_t>(const std::string& filename, MapMatrix& m);
};

/*Implementation of the functions of the Map Matrix class*/
template<typename value_t>
MapMatrix<value_t> &MapMatrix<value_t>::operator=(const MapMatrix<value_t> &m)
{
    nr = m.nr;
    nc = m.nc;
    data = m.data;

    return *this;
}

/*Insert a value into the map matrix*/
template<typename value_t>
void MapMatrix<value_t>::insert(const int &i, const int &j, const value_t &v)
{
    assert (i < nr && j < nc);
    data.insert_or_assign(std::make_tuple(i,j), v);
}

/*Output stream for the map matrix*/
template<typename value_t>
std::ostream &operator<<(std::ostream &o, MapMatrix<value_t> &m)
{
    for(int i = 0; i < m.nr; i++)
    {
        for(int j = 0; j < m.nc; j++)
        {
            auto val = m.data.find(std::make_tuple(i,j));
            if(val != m.data.end())
            {
                o << val->second << "\t";
            }
            else
            {
                o << 0.0 << "\t";
            }
        }
        o << std::endl;
    }
    return o;
}

/*Add operator to add two map mnatrices*/
template<typename value_t>
MapMatrix<value_t> operator+(const MapMatrix<value_t> &m1, const MapMatrix<value_t> &m2)
{
    assert(m1.nr == m2.nr && m1.nc == m2.nc);

    MapMatrix<value_t> res(m1);

    for(auto elem_m2: m2.data)
    {
        auto elem_m1 = res.data.find(elem_m2.first);
        if(elem_m1 == res.data.end())
        {
            res.data.insert(elem_m2);
        }
        else
        {
            res.data.insert_or_assign(elem_m1->first, elem_m1->second+elem_m2.second);
        }
    }
    return res;
}

template<typename value_t>
MapMatrix<value_t> &MapMatrix<value_t>::operator+=(const MapMatrix<value_t> &m)
{
    assert(this->nr == m.nr && this->nc == m.nc);

    for(auto elem_m: m.data)
    {
        auto elem_this = data.find(elem_m.first);
        if(elem_this == data.end())
        {
            data.insert(elem_m);
        }
        else
        {
            data.insert_or_assign(elem_this->first, elem_this->second+elem_m.second);
        }
    }
    return *this;
}

/*Multiplication of two map matrices. Using an initial transpose of the first matrix and the inherrrent soting of the underlying map type the multiplication is performed.*/
template<typename value_t>
MapMatrix<value_t> operator*(const MapMatrix<value_t> &m1, const MapMatrix<value_t> &m2)
{
    assert(m1.nc==m2.nr);
    MapMatrix<value_t> res(m1.nr, m2.nc);

    std::map<std::tuple<int, int>, value_t> transpose;
    for(const auto &[idx, val] : m1.data)
    {
        transpose.insert(std::make_pair(std::make_tuple(std::get<1>(idx), std::get<0>(idx)), val));
    }
 
    auto m2_iter = m2.data.begin();
        
    for (const auto &[idx, val] : transpose)
    {
        //Catch up with the first matrix
        while(std::get<0>(m2_iter->first) < std::get<1>(idx)) m2_iter++;

        //Let the first matrix catch up
        if(std::get<0>(m2_iter->first) > std::get<1>(idx)) continue;

        auto drag = m2_iter;
        while(std::get<0>(m2_iter->first) == std::get<0>(idx))
        {
            //Element present in both matrices
            auto res_idx = std::make_tuple(std::get<1>(idx), std::get<1>(m2_iter->first));
            auto current_elem = res.data.find(res_idx);
            if(current_elem == res.data.end())
            {
                res.data.emplace(res_idx, val*m2_iter->second);
            }
            else
            {
                res.data.insert_or_assign(res_idx, current_elem->second+(val*m2_iter->second));
            }
            m2_iter++;
        }
        m2_iter = drag;
    }

    return res;
}


/*Multiply the elements of the map matrix by a value a*/
template<typename value_t>
MapMatrix<value_t> operator*(const value_t &a, const MapMatrix<value_t> &m)
{
    MapMatrix<value_t> res(m);
    for (auto &elem : m.data)
    {
        res.data.insert_or_assign(elem.first, a*elem.second);
    }

    return res;
}

template<typename value_t>
std::vector<value_t> operator*(const MapMatrix<value_t> &m, const std::vector<value_t> &v)
{
    assert(v.size() == m.nc);

    std::vector<value_t> res(m.nr, 0.0);

    for (const auto &[idx, val] : m.data)
    {
        res[std::get<0>(idx)] +=  val*v[std::get<1>(idx)];
    }

    return res;
}

template<typename value_t>
MapMatrix<value_t> &MapMatrix<value_t>::operator*=(const MapMatrix<value_t> &m)
{
    assert(nc==m.nr);
    
    std::map<NxN, value_t> tmp;
    std::map<NxN, value_t> transpose;
    for(const auto &[idx, val] : data)
    {
        transpose.insert(std::make_pair(std::make_tuple(std::get<1>(idx), std::get<0>(idx)), val));
    }
 
    auto m_iter = m.data.begin();
        
    for (const auto &[idx, val] : transpose)
    {
        //Catch up with the first matrix
        while(std::get<0>(m_iter->first) < std::get<1>(idx)) m_iter++;

        //Let the first matrix catch up
        if(std::get<0>(m_iter->first) > std::get<1>(idx)) continue;

        auto drag = m_iter;
        while(std::get<0>(m_iter->first) == std::get<0>(idx))
        {
            //Element present in both matrices
            auto res_idx = std::make_tuple(std::get<1>(idx), std::get<1>(m_iter->first));
            auto current_elem = tmp.find(res_idx);
            if(current_elem == tmp.end())
            {
                tmp.emplace(res_idx, val*m_iter->second);
            }
            else
            {
                tmp.insert_or_assign(res_idx, current_elem->second+(val*m_iter->second));
            }
            m_iter++;
        }
        m_iter = drag;
    }
    nc = m.nc;
    data = tmp;

    return *this;
}


template<typename value_t>
MapMatrix<value_t> &MapMatrix<value_t>::operator*=(const double &a)
{
    for (auto &elem : data)
    {
        data.insert_or_assign(elem.first, a*elem.second);
    }
    return *this;
}

/*General template to load a map matrix from a file*/
template<typename value_t>
MapMatrix<value_t> LoadMapMatrix(const std::string& filename)
{
    int nr,nc;

    std::ifstream file(filename);

    /*Read the first line corresponding to the number of rows and columns*/
    file >> nr >> nc;
    
    MapMatrix<value_t> m(nr, nc);

    int j, k;
    value_t v;
    while(file >> j >> k >> v)
    {
        m.insert(j, k, v);
    }

    file.close();
    return m;
}

/*Specialization to load a complex map matrix from a file*/
template <>
inline MapMatrix<Cplx> LoadMapMatrix(const std::string& filename)
{
    int nr,nc;
    std::ifstream file(filename);

    file >> nr >> nc;

    MapMatrix<Cplx> m(nr,nc);

    int j,k;
    double v_real, v_imag;

    while(file >> j >> k >> v_real >> v_imag)
    {
        m.insert(j,k, Cplx(v_real, v_imag));
    }

    file.close();
    return m;
}

/*General template to write a map matrix to a file*/
template<typename value_t>
void Write(const std::string& filename, MapMatrix<value_t>& m)
{
    std::ofstream file(filename);

    file << m.nr << " " << m.nc << std::endl;

    for(const auto &elem : m.data)
    {
        file << std::get<0>(elem.first) << " " << std::get<1>(elem.first) << " " << elem.second << std::endl;
    }

    file.close();
}

/*Specialization to write complex valued map matrices to a file*/
template <>
inline void Write(const std::string &filename, MapMatrix<Cplx>& m)
{
    std::ofstream file(filename);

    file << m.nr << " " << m.nc << std::endl;

    for(const auto &elem : m.data)
    {
        file << std::get<0>(elem.first) << " " << std::get<1>(elem.first) << " " << std::real(elem.second) << " " << std::imag(elem.second) << std::endl;
    }

    file.close();
}