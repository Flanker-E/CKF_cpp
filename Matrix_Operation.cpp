#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include <iomanip>
#include "Matrix_Operation.h"
using namespace std;


//根据用户需求新建矩阵
Matrix InputMatrix()
{
	Matrix m;
	int row, col;
	cout << "输入行数与列数：" << endl;
	cin >> row >> col;
	double **enterMatrix;
	enterMatrix = (double **)malloc(row * sizeof(double *));
	for (int i = 0; i < row; i++)
		enterMatrix[i] = (double *)malloc(col * sizeof(double));
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
Matrix ZerosMatrix(int row, int col)
{
	Matrix m;
	double** matrix;
	matrix = (double**)malloc(row * sizeof(double*));
	for (int i = 0; i < row; i++)
		matrix[i] = (double*)malloc(col * sizeof(double));
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
Matrix TransposeMatrix(Matrix m, int free)
{
	int row = m.row;
	int col = m.col;
	Matrix mt = ZerosMatrix(col, row); //col和row反过来
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			mt.matrix[j][i] = m.matrix[i][j];
		}
	}
	if (free)
		freeMatrix(m);
	return mt;
}

//类似matlab的eye函数
Matrix EyeMatrix(int k)
{
	Matrix m = ZerosMatrix(k, k);
	for (int i = 0; i < k; i++)
		m.matrix[i][i] = 1;
	return m;
}

Matrix OnesMatrix(int row, int col, double num)
{
	Matrix m;
	double** matrix;
	matrix = (double**)malloc(row * sizeof(double*));
	for (int i = 0; i < row; i++)
		matrix[i] = (double*)malloc(col * sizeof(double));
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			matrix[i][j] = num;
		}
	}
	m.col = col;
	m.row = row;
	m.matrix = matrix;
	return m;
}

//新建矩阵并将相加值以新建矩阵返回,如果free=1（默认为0）则把输入矩阵的空间释放
Matrix add(Matrix m1, Matrix m2, int free)
{
	Matrix result = ZerosMatrix(m1.row, m1.col);
	for (int i = 0; i < m1.row; i++)
	{
		for (int j = 0; j < m1.col; j++)
		{
			result.matrix[i][j] = m1.matrix[i][j] + m2.matrix[i][j];
		}
	}
	if (free)
	{
		if ((free == 10) | (free == 11))
			freeMatrix(m1);
		if ((free == 01) | (free == 11))
			freeMatrix(m2);
	}
	return result;
}

//新建矩阵并将相减值以新建矩阵返回,如果free=1（默认为0）则把输入矩阵的空间释放
Matrix sub(Matrix m1, Matrix m2, int free)
{
	Matrix result = ZerosMatrix(m1.row, m1.col);
	for (int i = 0; i < m1.row; i++)
	{
		for (int j = 0; j < m1.col; j++)
		{
			result.matrix[i][j] = m1.matrix[i][j] - m2.matrix[i][j];
		}
	}
	if (free)
	{
		if ((free == 10) | (free == 11))
			freeMatrix(m1);
		if ((free == 01) | (free == 11))
			freeMatrix(m2);
	}
	return result;
}


//行列相乘
double calRowCol(Matrix M1, Matrix M2, int row, int col) //row为M1的行 col为m2的列
{
	double result = 0;
	int same = M1.col;
	for (int j = 0; j < same; j++)
	{
		result += M1.matrix[row][j] * M2.matrix[j][col];
	}

	return result;
}

//矩阵叉乘
Matrix Mul(Matrix m1, Matrix m2, int free)
{
	Matrix result = ZerosMatrix(m1.row, m2.col);
	for (int i = 0; i < m1.row; i++)
	{
		for (int j = 0; j < m2.col; j++)
		{
			result.matrix[i][j] = calRowCol(m1, m2, i, j);
		}
	}
	if (free)
	{
		if ((free == 10) | (free == 11))
			freeMatrix(m1);
		if ((free == 01) | (free == 11))
			freeMatrix(m2);
	}
	return result;
}

//矩阵数乘
Matrix numMul(Matrix m, double num, int free)
{
	Matrix result = ZerosMatrix(m.row, m.col);
	//cout << "数值:" << num << endl;
	for (int i = 0; i < m.row; i++)
	{
		for (int j = 0; j < m.col; j++)
		{
			result.matrix[i][j] = m.matrix[i][j] * num;
		}
	}
	if (free)
		freeMatrix(m);
	return result;
}

int isSymmetric(Matrix m)
{
	if (m.col == m.row)
	{
		for (int i = 0; i < m.col - 1; i++)
			for (int j = i + 1; j < m.col; j++)
			{
				if (m.matrix[i][j] != m.matrix[j][i])
					return false;
			}
		return true;
	}
	else
		return false;
}

//similar to matlab's diag input one row/col matrix and return diag.
Matrix diagMatrix(Matrix m, int free)
{
	int num = (m.row > m.col) ? m.row : m.col;
	Matrix diag = ZerosMatrix(num, num);
	for (int i = 0; i < num; i++)
		if (m.row > m.col)
			diag.matrix[i][i] = m.matrix[i][0];
		else
			diag.matrix[i][i] = m.matrix[0][i];
	if (free)
		freeMatrix(m);
	return diag;
}

//similar to matlab's remap
Matrix repmatMatrix(Matrix m, int rownum, int colnum, int free)
{
	Matrix result = ZerosMatrix((m.row*rownum), (m.col*colnum));
	for (int i = 0; i < rownum; i++)
		for (int j = 0; j < colnum; j++)
			passMatrix(result, m, i*m.row + 1, (i + 1)*m.row, j*m.col + 1, (j + 1)*m.col);
	if (free)
		freeMatrix(m);
	return result;
}

//Cholsky Decompose
Matrix Cholesky(Matrix ma)
{
	int k = ma.col;
	double sumi = 0;
	double sumj = 0;
	Matrix ml = ZerosMatrix(ma.row, ma.col);
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
	Matrix mY = ZerosMatrix(ma.row, ma.col);//存放Y结果的空矩阵
	Matrix mI = EyeMatrix(k);//创建一个单位阵
	/*L*Y=I*/
	double sumi = 0;
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
	Matrix mlt = ZerosMatrix(ma.row, ma.col);
	mlt = TransposeMatrix(ml);//求L的转置
	Matrix mX = ZerosMatrix(ma.row, ma.col);//存放X结果的空矩阵
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

//传递矩阵值给另一个矩阵
void passMatrix(Matrix to, Matrix from,
	int to_rowbegin, int to_rowend, int to_colbegin, int to_colend,
	int from_rowbegin, int from_rowend, int from_colbegin, int from_colend)
{
	////redirection
	if (to_rowend == 0)
		to_rowend = to.row;
	if (to_colend == 0)
		to_colend = to.col;
	if (from_rowend == 0)
		from_rowend = from.row;
	if (from_colend == 0)
		from_colend = from.col;

	for (int i = 0; i < (from_rowend - from_rowbegin + 1); i++)
		for (int j = 0; j < (from_colend - from_colbegin + 1); j++)
			to.matrix[to_rowbegin + i - 1][to_colbegin + j - 1] = from.matrix[from_rowbegin + i - 1][from_colbegin + j - 1];
}

//obtain mean value
double meanMatrix(Matrix m, int rowbegin, int rowend, int colbegin, int colend, int free)
{
	if (rowend == 0)
		rowend = m.row;
	if (colend == 0)
		colend = m.col;

	double mean = 0;
	int num = 0;
	for (int i = 0; i < (rowend - rowbegin + 1); i++)
		for (int j = 0; j < (colend - colbegin + 1); j++)
		{
			mean += m.matrix[rowbegin + i - 1][colbegin + j - 1];
			num++;
		}
	if (free)
		freeMatrix(m);
	mean = mean / num;
	return mean;
}

void sort(Matrix m, int row, int col)
{
	int nums[10] = { 4, 5, 2, 10, 7, 1, 8, 3, 6, 9 };
	int i, j, isSorted;
	double temp;
	if ((row != 0) && (col == 0))//sort row
	{
		//优化算法：最多进行 n-1 轮比较
		for (i = 0; i < m.col - 1; i++)
		{
			isSorted = 1;  //假设剩下的元素已经排序好了
			for (j = 0; j < m.col - 1 - i; j++) {
				if (m.matrix[row - 1][j] > m.matrix[row - 1][j + 1]) {
					temp = m.matrix[row - 1][j];
					m.matrix[row - 1][j] = m.matrix[row - 1][j + 1];
					m.matrix[row - 1][j + 1] = temp;
					isSorted = 0;  //一旦需要交换数组元素，就说明剩下的元素没有排序好
				}
			}
			if (isSorted) break; //如果没有发生交换，说明剩下的元素已经排序好了
		}
	}
	else//sort col
	{
		//优化算法：最多进行 n-1 轮比较
		for (i = 0; i < m.row - 1; i++)
		{
			isSorted = 1;  //假设剩下的元素已经排序好了
			for (j = 0; j < m.row - 1 - i; j++) {
				if (m.matrix[j][col - 1] > m.matrix[j + 1][col - 1]) {
					temp = m.matrix[j][col - 1];
					m.matrix[j][col - 1] = m.matrix[j + 1][col - 1];
					m.matrix[j + 1][col - 1] = temp;
					isSorted = 0;  //一旦需要交换数组元素，就说明剩下的元素没有排序好
				}
			}
			if (isSorted) break; //如果没有发生交换，说明剩下的元素已经排序好了
		}
	}
}

//return median value.row/col which input is 0 is null.
double medianMatrix(Matrix m_median, int row, int col)
{
	Matrix m = m_median;
	double median = 0;
	if ((row == 0) && (col == 0))
		return 0;
	else
	{
		sort(m, row, col);
		if ((row != 0) && (col == 0))//sort row
			if ((m.col % 2) == 0)
				median = (m.matrix[row - 1][m.col / 2 - 1] + m.matrix[row - 1][m.col / 2]) / 2.0;
			else
				median = m.matrix[row - 1][(int)floor((m.col) / 2.0)];
		else//sort col
			if ((m.row % 2) == 0)
				median = (m.matrix[m.row / 2 - 1][col - 1] + m.matrix[m.row / 2][col - 1]) / 2.0;
			else
				median = m.matrix[(int)floor((m.row) / 2.0)][col - 1];
	}
	//printMatrix(m);
	freeMatrix(m);
	return median;
}


//打印矩阵
void printMatrix(Matrix m)
{
	for (int i = 0; i < m.row; i++)
	{
		for (int j = 0; j < m.col; j++)
		{
			cout << setprecision(16) << m.matrix[i][j] << "  ";
		}
		cout  << endl;
	}
	cout << endl;
}

void freeMatrix(Matrix m)
{
	for (int n = 0; n < m.row; n++)
	{
		free(m.matrix[n]);
	}
	free(m.matrix);
}