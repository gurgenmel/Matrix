#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include "matrix.h"

#define EPS 1e-10
using namespace std;
using std::ostream;  using std::istream;  using std::endl;
using std::domain_error;



Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols)
{
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = 0;
        }
    }
}

Matrix::Matrix() : rows_(1), cols_(1)
{
    allocSpace();
    p[0][0] = 0;
}

Matrix::~Matrix()
{
    for (int i = 0; i < rows_; ++i) {
        delete[] p[i];
    }
    delete[] p;
}

Matrix::Matrix(const Matrix& m) : rows_(m.rows_), cols_(m.cols_)
{
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
}

Matrix& Matrix::operator=(const Matrix& m)
{
    if (this == &m) {
        return *this;
    }

    if (rows_ != m.rows_ || cols_ != m.cols_) {
        for (int i = 0; i < rows_; ++i) {
            delete[] p[i];
        }
        delete[] p;

        rows_ = m.rows_;
        cols_ = m.cols_;
        allocSpace();
    }

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator+=(const Matrix& m)
{
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			p[i][j] += m.p[i][j];
		}
	}
	return *this;
}

Matrix& Matrix::operator-=(const Matrix& m)
{
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			p[i][j] -= m.p[i][j];
		}
	}
	return *this;
}

Matrix& Matrix::operator*=(const Matrix& m)
{
	Matrix x;
	if (cols_ == m.rows_)
	{
		Matrix temp(rows_, m.cols_);
		for (int i = 0; i < temp.rows_; ++i) {
			for (int j = 0; j < temp.cols_; ++j) {
				for (int k = 0; k < cols_; ++k) {
					temp.p[i][j] += (p[i][k] * m.p[k][j]);
				}
			}
		}
		return (*this = temp);
	}
	else
	{
		return x;
	}
}

Matrix& Matrix::operator*=(double num)
{
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			p[i][j] *= num;
		}
	}
	return *this;
}





Matrix Matrix::AdProduct(const Matrix& m)
{
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			p[i][j] *= m.p[i][j];
		}
	}
	return *this;
}
Matrix Matrix::getCofactor(Matrix& temp , int f, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != f && col != q)
			{
				temp.p[i][j++] = p[row][col];

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
	return temp;
}
double Matrix::det(const Matrix& m,int n)
{
	int D = 0; // Initialize result 
			   //  Base case : if matrix contains single element 
	if (n == 1)
		return m.p[0][0];

	Matrix temp(n-1,n-1); // To store cofactors 

	int sign = 1;  // To store sign multiplier 

				   // Iterate for each element of first row 
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of mat[0][f] 
		temp=getCofactor(temp, 0, f, n);
		D += sign * p[0][f] * det(temp, n - 1);

		// terms are to be added with alternate sign 
		sign = -sign;
	}

	return D;
}

void Matrix::allocSpace()
{
    p = new double*[rows_];
    for (int i = 0; i < rows_; ++i) {
        p[i] = new double[cols_];
    }
}





Matrix operator+(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp += m2);
}

Matrix operator-(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp -= m2);
}

Matrix operator*(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp *= m2);
}

Matrix operator*(const Matrix& m, double num)
{
    Matrix temp(m);
    return (temp *= num);
}

Matrix operator*(double num, const Matrix& m)
{
    return (m * num);
}


ostream& operator<<(ostream& os, const Matrix& m)
{
    for (int i = 0; i < m.rows_; ++i) {
        os << m.p[i][0];
        for (int j = 1; j < m.cols_; ++j) {
            os << " " << m.p[i][j];
        }
        os << endl;
    }
    return os;
}

istream& operator>>(istream& is, Matrix& m)
{
    for (int i = 0; i < m.rows_; ++i) {
        for (int j = 0; j < m.cols_; ++j) {
            is >> m.p[i][j];
        }
    }
    return is;
}
double Matrix::dotProduct(Matrix a, Matrix b)
{
	double sum = 0;
	for (int i = 0; i < a.rows_; ++i) {
		sum += (a(i, 0) * b(i, 0));
	}
	return sum;
}
double Matrix::trace(Matrix a)
{
	double sum = 0;
	for (int i = 0; i < a.rows_; ++i) {
		for (int j = 0; j < a.cols_; j++)
			if (i == j)
				sum += a(i, j);
	}
	return sum;
}
double Matrix::exterior(Matrix a,Matrix b)
{
	double sum = 0;
	sum = (a.dotProduct(a, b) - a.dotProduct(b, a))/2;
	return sum;
}
double Matrix::mnorm(Matrix a)
{
	double max = 0;
	for (int i = 0; i < a.rows_; ++i) {
		if (a.p[i][0] > max)
			max = a.p[i][0];
	}
	return max;
}
double Matrix::enorm(Matrix a)
{
	double max = 0;
	for (int i = 0; i < a.rows_; ++i) {
		for(int j=0;j<a.cols_;j++)
			max += a.p[i][j]*a.p[i][j];
	}
	max = sqrt(max);
	return max;
}
double Matrix::angel(Matrix a, Matrix b)
{
	double ang = 0;
	ang = (a.dotProduct(a, b)) / (a.enorm(a)*b.enorm(b));
	ang = acos(ang);
	return ang;
}
int Matrix::rank(Matrix a)
{
	int rank = max(a.rows_, a.cols_);
	vector<char> line_used(a.rows_);
	for (int i = 0; i<a.cols_; ++i) {
		int j;
		for (j = 0; j<a.rows_; ++j)
			if (!line_used[j] && abs(a.p[j][i]) > EPS)
				break;
		if (j == a.rows_)
			--rank;
		else {
			line_used[j] = true;
			for (int l = i + 1; l<a.cols_; l++)
				a.p[j][l] /= a.p[j][i];
			for (int k = 0; k<a.rows_; k++)
				if (k != j && abs(a.p[k][i]) > EPS)
					for (int l = i + 1; l<a.cols_; l++)
						a.p[k][l] -= a.p[j][l] * a.p[k][i];
		}
	}
	return rank;
}
Matrix Matrix::transpose()
{
	Matrix ret(cols_, rows_);
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			ret.p[j][i] = p[i][j];
		}
	}
	return ret;
}
Matrix Matrix::inverse()
{
	Identity I(rows_,rows_);
	Matrix AI = Matrix::augment(*this, I);
	Matrix U = AI.gaussianEliminate();
	Matrix IAInverse = U.rowReduceFromGaussian();
	Matrix AInverse(rows_, cols_);
	for (int i = 0; i < AInverse.rows_; ++i) {
		for (int j = 0; j < AInverse.cols_; ++j) {
			AInverse(i, j) = IAInverse(i, j + cols_);
		}
	}
	return AInverse;
}
Matrix Matrix::augment(Matrix A, Matrix B)
{
	Matrix AB(A.rows_, A.cols_ + B.cols_);
	for (int i = 0; i < AB.rows_; ++i) {
		for (int j = 0; j < AB.cols_; ++j) {
			if (j < A.cols_)
				AB(i, j) = A(i, j);
			else
				AB(i, j) = B(i, j - B.cols_);
		}
	}
	return AB;
}

Matrix Matrix::gaussianEliminate()
{
	Matrix Ab(*this);
	int rows = Ab.rows_;
	int cols = Ab.cols_;
	int Acols = cols - 1;

	int i = 0; // row tracker
	int j = 0; // column tracker

			   // iterate through the rows
	while (i < rows)
	{
		// find a pivot for the row
		bool pivot_found = false;
		while (j < Acols && !pivot_found)
		{
			if (Ab(i, j) != 0) { // pivot not equal to 0
				pivot_found = true;
			}
			else { // check for a possible swap
				int max_row = i;
				double max_val = 0;
				for (int k = i + 1; k < rows; ++k)
				{
					double cur_abs = Ab(k, j) >= 0 ? Ab(k, j) : -1 * Ab(k, j);
					if (cur_abs > max_val)
					{
						max_row = k;
						max_val = cur_abs;
					}
				}
				if (max_row != i) {
					Ab.swapRows(max_row, i);
					pivot_found = true;
				}
				else {
					j++;
				}
			}
		}

		// perform elimination as normal if pivot was found
		if (pivot_found)
		{
			for (int t = i + 1; t < rows; ++t) {
				for (int s = j + 1; s < cols; ++s) {
					Ab(t, s) = Ab(t, s) - Ab(i, s) * (Ab(t, j) / Ab(i, j));
					if (Ab(t, s) < EPS && Ab(t, s) > -1 * EPS)
						Ab(t, s) = 0;
				}
				Ab(t, j) = 0;
			}
		}

		i++;
		j++;
	}

	return Ab;
}

Matrix Matrix::rowReduceFromGaussian()
{
	Matrix R(*this);
	int rows = R.rows_;
	int cols = R.cols_;

	int i = rows - 1; // row tracker
	int j = cols - 2; // column tracker

					  // iterate through every row
	while (i >= 0)
	{
		// find the pivot column
		int k = j - 1;
		while (k >= 0) {
			if (R(i, k) != 0)
				j = k;
			k--;
		}

		// zero out elements above pivots if pivot not 0
		if (R(i, j) != 0) {

			for (int t = i - 1; t >= 0; --t) {
				for (int s = 0; s < cols; ++s) {
					if (s != j) {
						R(t, s) = R(t, s) - R(i, s) * (R(t, j) / R(i, j));
						if (R(t, s) < EPS && R(t, s) > -1 * EPS)
							R(t, s) = 0;
					}
				}
				R(t, j) = 0;
			}

			// divide row by pivot
			for (int k = j + 1; k < cols; ++k) {
				R(i, k) = R(i, k) / R(i, j);
				if (R(i, k) < EPS && R(i, k) > -1 * EPS)
					R(i, k) = 0;
			}
			R(i, j) = 1;

		}

		i--;
		j--;
	}

	return R;
}
void Matrix::swapRows(int r1, int r2)
{
	double *temp = p[r1];
	p[r1] = p[r2];
	p[r2] = temp;
}
