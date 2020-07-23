// Matrix_test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include "Matrix_Operation.h"
using namespace std;


int main()
{
    int num = 0;
    do
    {
        cout << "*************************************\n";
        cout << "*              菜单                 *\n";
        cout << "*          0.cholesky               *\n";
        cout << "*          1.矩阵相加               *\n";
        cout << "*          2.矩阵相减               *\n";
        cout << "*          3.矩阵相乘               *\n";
        cout << "*          4.矩阵数乘               *\n";
        cout << "*          5.test_free              *\n";
        cout << "*************************************\n";
        cin >> num;
        if (1 == num || 2 == num || 3 == num)
        {
            cout << "请输入矩阵1" << endl;
            Matrix m1 = CreateMatrix();
            cout << "请输入矩阵2" << endl;
            Matrix m2 = CreateMatrix();
            cout << "两矩阵为" << endl;
            printMatrix(m1);
            cout << endl;
            printMatrix(m2);
            /*cout << m1.matrix;*/
            switch (num)
            {
            case 1:
            {
                if (m1.col != m2.col || m1.row != m2.row)
                {
                    cout << "行列不同" << endl;
                }
                else
                {
                    cout << "结果为：" << endl;
                    printMatrix(add(m1, m2));
                }
                break;
            }
            case 2:
            {

                if (m1.col != m2.col || m1.row != m2.row)
                {
                    cout << "参数错误" << endl;
                }
                else
                {
                    cout << "结果为：" << endl;
                    printMatrix(sub(m1, m2));
                }
                break;
            }
            case 3:
            {
                if (m1.col != m2.row)
                {
                    cout << "参数错误" << endl;
                }
                else
                {
                    cout << "结果为：" << endl;
                    printMatrix(Mul(m1, m2));
                }
                break;
            }
            default:
                break;
            }
            free(m1.matrix);
            free(m2.matrix);
        }
        else if (4 == num)
        {
            int number = 1;
            cout << "请输入矩阵" << endl;
            Matrix m = CreateMatrix();
            cout << "请输入数值" << endl;
            cin >> number;
            cout << "矩阵为：" << endl;
            printMatrix(m);
            cout << "数值为：" << endl;
            cout << number << endl;
            printMatrix(numMul(m, number));
        }
        else if (5 == num)
        {
            for (int n = 0; n < 1000; n++)
            {
                int col = 5; int row = 5;
                Matrix m1 = InitMatrix(row, col);
                freeMatrix(m1);
            } 
        }
        else if (0 == num) //cholesky
        {
            ////////Cholesky////////////
            Matrix ma = CreateMatrix();
            if (ma.col == ma.row)
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

                printMatrix(ml);
                ////////////inverse///////////////
                Matrix mY = InitMatrix(ma.row, ma.col);
                Matrix mI = EyeMatrix(k);
                /*L*Y=I*/
                for (int j = 0; j < k; j++)
                {
                    for (int i = 0; i < k; i++)
                    {
                        if (i > 0)
                        {
                            for (int n = 0; n < i ; n++)
                                sumi += ml.matrix[i][n] * mY.matrix[n][j];
                        }
                        mY.matrix[i][j] = (mI.matrix[i][j] - sumi) / ml.matrix[i][i];
                        sumi = 0;
                    }
                }
                /*LT*X=Y*/
                Matrix mlt = InitMatrix(ma.row, ma.col);
                mlt = TransposeMatrix(ml, ml.row, ml.col);
                Matrix mX = InitMatrix(ma.row, ma.col);
                for (int j = k - 1; j >= 0; j--)
                {
                    for (int i = k - 1; i >= 0; i--)
                    {
                        if (i < k - 1)
                        {
                            for (int n = i + 1; n < k ;n++)
                                sumi += mlt.matrix[i][n] * mX.matrix[n][j];
                        }
                        mX.matrix[i][j] = (mY.matrix[i][j] - sumi) / mlt.matrix[i][i];
                        sumi = 0;
                    }
                }
                printMatrix(mX);
            }
        }
        cout << "按回车继续....";
        /*cout << m1.matrix;*/
        /*free(m1.matrix);
        free(m2);*/
        getchar();
        getchar();
        system("cls");
    } while (0 == num || 1 == num || 2 == num || 3 == num || 4 == num || 5 == num);
    return 0;
}

//int main()
//{
//    std::cout << "Hello World!\n";
//}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started:
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
