#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include "Matrix_Operation.h"
using namespace std;


//�����û������½�����
Matrix CreateMatrix()
{
	Matrix m;
	int row, col;
	cout << "����������������" << endl;
	cin >> row >> col;
	float **enterMatrix;
	enterMatrix = (float **)malloc(row * sizeof(float *));
	for (int i = 0; i < row; i++)
		enterMatrix[i] = (float *)malloc(col * sizeof(float));
	cout << "������ľ���" << endl;
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

//��ʼ��һ����Ϊrow��Ϊcol����
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

//ȡת�ã�������һ����matrix
Matrix TransposeMatrix(Matrix m, int row, int col)
{
	Matrix mt = InitMatrix(col, row); //col��row������
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			mt.matrix[j][i] = m.matrix[i][j];
		}
	}
	return mt;
}

//����matlab��eye����
Matrix EyeMatrix(int k)
{
	Matrix m = InitMatrix(k, k);
	for (int i = 0; i < k; i++)
		m.matrix[i][i] = 1;
	return m;
}


//�½����󲢽����ֵ���½����󷵻�
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

//�½����󲢽����ֵ���½����󷵻�
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

//�������
float calRowCol(Matrix M1, Matrix M2, int row, int col) //rowΪM1���� colΪm2����
{
	float result = 0;
	int same = M1.col;
	for (int j = 0; j < same; j++)
	{
		result += M1.matrix[row][j] * M2.matrix[j][col];
	}

	return result;
}

//������
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

//��������
Matrix numMul(Matrix m, int num)
{
	Matrix result = InitMatrix(m.row, m.col);
	cout << "��ֵ:" << num << endl;
	for (int i = 0; i < m.row; i++)
	{
		for (int j = 0; j < m.col; j++)
		{
			result.matrix[i][j] = m.matrix[i][j] * num;
		}
	}
	return result;
}

//��ӡ����
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