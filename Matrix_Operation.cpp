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

bool isSymmetric(Matrix m)
{
	if (m.col == m.row)
	{
		for (int i = 0; i < m.col-1; i++)
			for (int j = i + 1; j < m.col; j++)
			{
				if (m.matrix[i][j]!= m.matrix[j][i])
					return false;
			}
		return true;
	}
	else
		return false;
}

//Cholsky Decompose
Matrix Cholesky(Matrix ma)
{
	int k = ma.col;
	float sumi = 0;
	float sumj = 0;
	Matrix ml = InitMatrix(ma.row, ma.col);
	for (int j = 0; j < k; j++)
	{
		if (j - 1 >= 0)
		{
			for (int n = 0; n < j; n++)
			{
				sumj += ml.matrix[j][n] * ml.matrix[j][n];
			}
		}
		ml.matrix[j][j] = sqrt(ma.matrix[j][j] - sumj); //对角线上的元素
		sumj = 0;
		for (int i = j + 1; i < k; i++)
		{
			if (j - 1 >= 0)
			{
				for (int n = 0; n < j; n++)
				{
					sumi += ml.matrix[j][n] * ml.matrix[i][n];
				}
			}
			ml.matrix[i][j] = (ma.matrix[i][j] - sumi) / (ml.matrix[j][j]); //其余同列元素
			sumi = 0;
		}
	}
	//printMatrix(ml);
	return ml;
}

//使用cholesky的求逆
Matrix CholeskyInverse(Matrix ma)
{
	////////////Cholesky////////////
	Matrix ml = Cholesky(ma);
	int k = ml.col;
	////////////inverse///////////////
	Matrix mY = InitMatrix(ma.row, ma.col);//存放Y结果的空矩阵
	Matrix mI = EyeMatrix(k);//创建一个单位阵
	/*L*Y=I*/
	float sumi = 0;
	for (int j = 0; j < k; j++)
	{
		for (int i = 0; i < k; i++)
		{
			if (i > 0)
			{
				for (int n = 0; n < i; n++)
					sumi += ml.matrix[i][n] * mY.matrix[n][j];
			}
			mY.matrix[i][j] = (mI.matrix[i][j] - sumi) / ml.matrix[i][i];
			sumi = 0;
		}
	}
	/*LT*X=Y*/
	Matrix mlt = InitMatrix(ma.row, ma.col);
	mlt = TransposeMatrix(ml, ml.row, ml.col);//求L的转置
	Matrix mX = InitMatrix(ma.row, ma.col);//存放X结果的空矩阵
	for (int j = k - 1; j >= 0; j--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			if (i < k - 1)
			{
				for (int n = i + 1; n < k; n++)
					sumi += mlt.matrix[i][n] * mX.matrix[n][j];
			}
			mX.matrix[i][j] = (mY.matrix[i][j] - sumi) / mlt.matrix[i][i];
			sumi = 0;
		}
	}
	//释放过程中创建的矩阵
	freeMatrix(ml);
	freeMatrix(mlt);
	freeMatrix(mY);
	freeMatrix(mI);
	//printMatrix(mX);
	return mX;
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