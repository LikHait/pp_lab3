#include "matrix.h"
#include <random>
#include "mpi.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dist(0, 100);
std::uniform_real_distribution<> dast(-10, 10);


std::ostream & operator<<(std::ostream & os, Elem & el)
{
	if (el.Re != 0)
	{
		os << el.Re;
		if (el.Im > 0)
			os << " + " << el.Im << 'i';
		if (el.Im < 0)
			os << " - " << (-1) * el.Im << 'i';
	}
	else
		if (el.Im != 0)
			os << el.Im << 'i';
	return os;
}

Elem Elem::operator*(const Elem & el)
{
	Elem tmp;
	tmp.Re = this->Re * el.Re - this->Im * el.Im;
	tmp.Im = this->Re * el.Im + this->Im * el.Re;
	return tmp;
}

Elem Elem::operator+(const Elem & el)
{
	return Elem(*this) += el;
}

Elem & Elem::operator+=(const Elem & el)
{
	this->Re += el.Re;
	this->Im += el.Im;
	return *this;
}

bool Elem::operator==(const Elem & el) const
{
	return (this->Re == el.Re) && (this->Im == el.Im);
}

bool Elem::operator!=(const Elem & el) const
{
	return !(*this == el);
}

//end Elem


void SparseMatrix::genMatrix(int den)
{
	for (int i = 0; i < _row; ++i)
		for (int j = 0; j < _col; ++j)
		{
			if (dist(gen) < den)
			{
				Elem el(dast(gen), 0); // dast(gen));
				v_elem.emplace_back(el);
				v_row.emplace_back(i);
				v_col.emplace_back(j);
			}
		}
}

SparseMatrix SparseMatrix::getTransposed()
{
	SparseMatrix tmp(_col, _row);
	for (int i = 0; i < _col; ++i)
		for (int j = 0; j < v_col.size(); ++j)
			if (v_col[j] == i)
			{
				tmp.v_elem.emplace_back(v_elem[j]);
				tmp.v_col.emplace_back(v_row[j]);
				tmp.v_row.emplace_back(v_col[j]);
			}
	return tmp;
}

SparseMatrix SparseMatrix::getTransposedV2() //для больших размеров матриц
{
	SparseMatrix tmp(_col, _row);
	std::vector<std::vector<int>> v1(_col);
	std::vector<std::vector<Elem>> v2(_col);
	int size = v_elem.size();
	for (int i = 0; i < size; ++i)
	{
		v1[v_col[i]].emplace_back(v_row[i]);
		v2[v_col[i]].emplace_back(v_elem[i]);
	}
	tmp.v_row.resize(size);
	int tmpSize = 0;
	for (int i = 0; i < _col; ++i)
	{
		tmp.v_col.insert(tmp.v_col.end(), v1[i].begin(), v1[i].end());
		tmp.v_elem.insert(tmp.v_elem.end(), v2[i].begin(), v2[i].end());
		std::fill(tmp.v_row.begin() + tmpSize, tmp.v_row.begin() + tmpSize + v1[i].size(), i);
		tmpSize += v1[i].size();
	}
	return tmp;
}

bool SparseMatrix::operator==(const SparseMatrix & sm) const
{
	if (this->v_elem.size() != sm.v_elem.size())
		return false;
	int size = this->v_elem.size();
	for (int i = 0; i < size; ++i)
		if (this->v_elem[i] != sm.v_elem[i] || this->v_col[i] != sm.v_col[i] || this->v_row[i] != sm.v_row[i])
			return false;
	return true;
}

std::ostream & operator<<(std::ostream & os, SparseMatrix & mx)
{
	for (int i = 0; i < mx.v_elem.size(); ++i)
		os << '(' << mx.v_row[i] << ',' << mx.v_col[i] << ") = " << mx.v_elem[i] << std::endl;
	return os;
}