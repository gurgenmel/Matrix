#pragma once



#include <iostream>

class Matrix {
public:
	Matrix(int, int);
	Matrix();
	~Matrix();
	Matrix(const Matrix&);
	Matrix& operator=(const Matrix&);
	Matrix& operator+=(const Matrix&);
	Matrix& operator-=(const Matrix&);
	Matrix& operator*=(const Matrix&);
	Matrix& operator*=(double);

	inline double& operator()(int x, int y) { return p[x][y]; }
	Matrix AdProduct(const Matrix& m);
	Matrix getCofactor( Matrix& temp, int p, int q, int n);
	double det(const Matrix&m, int n);
	static double dotProduct(Matrix, Matrix);
	static double trace(Matrix);
	static double exterior(Matrix, Matrix);
	static double mnorm(Matrix);
	static double enorm(Matrix);
	static double angel(Matrix, Matrix);
	static int rank(Matrix);
	friend std::ostream& operator<<(std::ostream&, const Matrix&);
	friend std::istream& operator>>(std::istream&, Matrix&);
	void swapRows(int, int);
	Matrix transpose();
	static Matrix augment(Matrix, Matrix);
	Matrix gaussianEliminate();
	Matrix rowReduceFromGaussian();
	Matrix inverse();

protected:
	int rows_, cols_;
	double **p;

	void allocSpace();
	
};

Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
class Identity : public Matrix
{
public:
	Identity(int rows,int cols) 
	{
		rows_ = rows;
		cols_=cols;
		allocSpace();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < rows; ++j) {
				if (i == j) {
					p[i][j] = 1;
				}
				else {
					p[i][j] = 0;
				}
			}
		}
	}
};
class Diagonal : public Matrix
{
public:
	Diagonal(int rows, int d)
	{
		rows_ = rows;
		cols_ = rows;
		allocSpace();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < rows; ++j) {
				if (i == j) {
					p[i][j] = d;
				}
				else {
					p[i][j] = 0;
				}
			}
		}
	}
};
class TriangleUp : public Matrix
{
public:
	TriangleUp(int n, int m)
	{
		rows_ = n;
		cols_ = m;
		allocSpace();
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				p[i][j] = 0;
			}
		}
	}
	friend std::istream& operator>>(std::istream& is, TriangleUp& m)
	{
		for (int i = 0; i < m.rows_; ++i) {
			for (int j = 0; j < m.cols_; ++j) {
				if(j>=i)
					is >> m.p[i][j];
			}
		}
		return is;
	}
};
class TriangleDown : public Matrix
{
public:
	TriangleDown(int n, int m)
	{
		rows_ = n;
		cols_ = m;
		allocSpace();
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				p[i][j] = 0;
			}
		}
	}
	friend std::istream& operator>>(std::istream& is, TriangleDown& m)
	{
		for (int i = 0; i < m.rows_; ++i) {
			for (int j = 0; j < m.cols_; ++j) {
				if (j<=i)
					is >> m.p[i][j];
			}
		}
		return is;
	}
};
class Sym : public Matrix
{
public:
	Sym(int n, int m)
	{
		rows_ = n;
		cols_ = m;
		allocSpace();
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				p[i][j] = 0;
			}
		}
	}
	friend std::istream& operator>>(std::istream& is, Sym& m)
	{
		for (int i = 0; i < m.rows_; ++i) {
			for (int j = 0; j < m.cols_; ++j) {
				if (j <= i)
					is >> m.p[i][j];
				//else
					//m.p[i][j] = p[j][i];
			}
		}
		
		return is;
	}
};

