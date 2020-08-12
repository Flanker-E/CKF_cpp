#pragma once

typedef struct
{
	//结构体
	int row, col;
	//二维指针，目的是动态分配内存
	double **matrix;
} Matrix;

Matrix InputMatrix();

Matrix ZerosMatrix(int row, int col);

Matrix TransposeMatrix(Matrix m, int free = false);

Matrix EyeMatrix(int k);

Matrix OnesMatrix(int row, int col, double num = 1);

Matrix add(Matrix m1, Matrix m2, int free = false);

Matrix sub(Matrix m1, Matrix m2, int free = false);

double calRowCol(Matrix M1, Matrix M2, int row, int col);

Matrix Mul(Matrix m1, Matrix m2, int free = false);

Matrix numMul(Matrix m, double num, int free = false);

int isSymmetric(Matrix m);

Matrix diagMatrix(Matrix m, int free = false);

Matrix repmatMatrix(Matrix m, int rownum, int colnum, int free = false);

Matrix Cholesky(Matrix ma);

Matrix CholeskyInverse(Matrix ma);

void passMatrix(Matrix to, Matrix from,
	int to_rowbegin = 1, int to_rowend = 0, int to_colbegin = 1, int to_colend = 0,
	int from_rowbegin = 1, int from_rowend = 0, int from_colbegin = 1, int from_colend = 0);

double meanMatrix(Matrix m, int rowbegin = 1, int rowend = 0, int colbegin = 1, int colend = 0, int free = false);

double medianMatrix(Matrix m, int row = 0, int col = 0);

void printMatrix(Matrix m);

void freeMatrix(Matrix m);