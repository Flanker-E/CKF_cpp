#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include "Matrix_Operation.h"
using namespace std;


//根据用户需求新建矩阵
Matrix CreateMatrix()
{
	Matrix m;
	int row, col;
	cout << "输入行数与列数：" << endl;
	cin >> row >> col;
	float **enterMatrix;
	enterMatrix = (float **)malloc(row * sizeof(float *));
	for (int i = 0; i < row; i++)
		enterMatrix[i] = (float *)malloc(col * sizeof(float));
	cout << "输入你的矩阵：" << endl;
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			cin >> enterMatrix[i][j];
		}
	}
	m.col = col;
	m.row = row;
	m.matrix = enterMatrix;
	return m;
}

//初始化一个行为row列为col矩阵
Matrix InitMatrix(int row, int col)
{
	Matrix m;
	float** matrix;
	matrix = (float**)malloc(row * sizeof(float*));
	for (int i = 0; i < row; i++)
		matrix[i] = (float*)malloc(col * sizeof(float));
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			matrix[i][j] = 0;
		}
	}
	m.col = col;
	m.row = row;
	m.matrix = matrix;
	return m;
}

//取转置，并返回一个新matrix
Matrix TransposeMatrix(Matrix m, int row, int col)
{
	Matrix mt = InitMatrix(col, row); //col和row反过来
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			mt.matrix[j][i] = m.matrix[i][j];
		}
	}
	return mt;
}

//类似matlab的eye函数
Matrix EyeMatrix(int k)
{
	Matrix m = InitMatrix(k, k);
	for (int i = 0; i < k; i++)
		m.matrix[i][i] = 1;
	return m;
}


//新建矩阵并将相加值以新建矩阵返回
Matrix add(Matrix m1, Matrix m2)
{
	Matrix result = InitMatrix(m1.row, m1.col);
	for (int i = 0; i < m1.row; i++)
	{
		for (int j = 0; j < m1.col; j++)
		{
			result.matrix[i][j] = m1.matrix[i][j] + m2.matrix[i][j];
		}
	}
	return result;
}

//新建矩阵并将相减值以新建矩阵返回
Matrix sub(Matrix m1, Matrix m2)
{
	Matrix result = InitMatrix(m1.row, m1.col);
	for (int i = 0; i < m1.row; i++)
	{
		for (int j = 0; j < m1.col; j++)
		{
			result.matrix[i][j] = m1.matrix[i][j] - m2.matrix[i][j];
		}
	}
	return result;
}

//行列相乘
float calRowCol(Matrix M1, Matrix M2, int row, int col) //row为M1的行 col为m2的列
{
	float result = 0;
	int same = M1.col;
	for (int j = 0; j < same; j++)
	{
		result += M1.matrix[row][j] * M2.matrix[j][col];
	}

	return result;
}

//矩阵叉乘
Matrix Mul(Matrix m1, Matrix m2)
{
	Matrix result = InitMatrix(m1.row, m2.col);
	for (int i = 0; i < m1.row; i++)
	{
		for (int j = 0; j < m2.col; j++)
		{
			result.matrix[i][j] = calRowCol(m1, m2, i, j);
		}
	}
	return result;
}

//矩阵数乘
Matrix numMul(Matrix m, int num)
{
	Matrix result = InitMatrix(m.row, m.col);
	cout << "数值:" << num << endl;
	for (int i = 0; i < m.row; i++)
	{
		for (int j = 0; j < m.col; j++)
		{
			result.matrix[i][j] = m.matrix[i][j] * num;
		}
	}
	return result;
}

//打印矩阵
void printMatrix(Matrix m)
{
	for (int i = 0; i < m.row; i++)
	{
		for (int j = 0; j < m.col; j++)
		{
			cout << m.matrix[i][j] << "  ";
		}
		cout << endl;
	}
}

void freeMatrix(Matrix m)
{
	for (int n = 0; n < m.row; n++)
	{
		free(m.matrix[n]);
	}
	free(m.matrix);
}