#pragma once
#include <vector>
#include <iostream>

struct Elem
{
	Elem(double re = 0, double im = 0): Re(re), Im(im) {}
	Elem operator*(const Elem &el);
	Elem operator+(const Elem &el);
	Elem& operator+=(const Elem &el);
	bool operator==(const Elem &el) const;
	bool operator!=(const Elem &el) const;

	friend std::ostream & operator<<(std::ostream & os, Elem & el);

	double Re;
	double Im;
};

class SparseMatrix
{
public:
	SparseMatrix(int row, int col) : _col(col), _row(row) {}
	SparseMatrix(const SparseMatrix &sm) : v_elem(sm.v_elem), v_row(sm.v_row), v_col(sm.v_col), _row(sm._row), _col(sm._col) {}
	void genMatrix(int den);
	SparseMatrix getTransposed();
	SparseMatrix getTransposedV2();

	void sortByRow();
	bool operator==(const SparseMatrix &sm) const;
	friend std::ostream & operator<<(std::ostream & os, SparseMatrix & mx);
public:
	std::vector<Elem> v_elem;
	std::vector<int> v_row;
	std::vector<int> v_col;
	int _row;
	int _col;
};

