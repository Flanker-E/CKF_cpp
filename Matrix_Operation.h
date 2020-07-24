#pragma once

typedef struct
{
	//结构体
	int row, col;
	//二维指针，目的是动态分配内存
	float **matrix;
} Matrix;

Matrix CreateMatrix();

Matrix InitMatrix(int row, int col);

Matrix TransposeMatrix(Matrix m, int row, int col);

Matrix EyeMatrix(int k);

Matrix add(Matrix m1, Matrix m2);

Matrix sub(Matrix m1, Matrix m2);

float calRowCol(Matrix M1, Matrix M2, int row, int col);

Matrix Mul(Matrix m1, Matrix m2);

Matrix numMul(Matrix m, int num);

bool isSymmetric(Matrix m);

Matrix Cholesky(Matrix ma);

Matrix CholeskyInverse(Matrix ma);

void printMatrix(Matrix m);

void freeMatrix(Matrix m);