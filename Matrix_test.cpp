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
        cout << "*              �˵�                 *\n";
        cout << "*          0.cholesky_inverse      *\n";
        //cout << "*          1.�������               *\n";
        //cout << "*          2.�������               *\n";
        //cout << "*          3.�������               *\n";
        //cout << "*          4.��������               *\n";
        cout << "*          1.test_free              *\n";
        cout << "*************************************\n";
        cin >> num;
		if (1 == num)
		{
			for (int n = 0; n < 1000; n++)//�ظ�1000�ο�free�ܲ�������ڴ�ռ��
			{
				int col = 5; int row = 5;
				Matrix m1 = InitMatrix(row, col);
				freeMatrix(m1);
			}
		}
		else if (0 == num) //cholesky
		{
			Matrix ma = CreateMatrix();

			if (isSymmetric(ma))
				Matrix mX = CholeskyInverse(ma);
			printMatrix(ma);
		}
        //if (1 == num || 2 == num || 3 == num)
        //{
        //    cout << "���������1" << endl;
        //    Matrix m1 = CreateMatrix();
        //    cout << "���������2" << endl;
        //    Matrix m2 = CreateMatrix();
        //    cout << "������Ϊ" << endl;
        //    printMatrix(m1);
        //    cout << endl;
        //    printMatrix(m2);
        //    /*cout << m1.matrix;*/
        //    switch (num)
        //    {
        //    case 1:
        //    {
        //        if (m1.col != m2.col || m1.row != m2.row)
        //        {
        //            cout << "���в�ͬ" << endl;
        //        }
        //        else
        //        {
        //            cout << "���Ϊ��" << endl;
        //            printMatrix(add(m1, m2));
        //        }
        //        break;
        //    }
        //    case 2:
        //    {
        //        if (m1.col != m2.col || m1.row != m2.row)
        //        {
        //            cout << "��������" << endl;
        //        }
        //        else
        //        {
        //            cout << "���Ϊ��" << endl;
        //            printMatrix(sub(m1, m2));
        //        }
        //        break;
        //    }
        //    case 3:
        //    {
        //        if (m1.col != m2.row)
        //        {
        //            cout << "��������" << endl;
        //        }
        //        else
        //        {
        //            cout << "���Ϊ��" << endl;
        //            printMatrix(Mul(m1, m2));
        //        }
        //        break;
        //    }
        //    default:
        //        break;
        //    }
        //    free(m1.matrix);
        //    free(m2.matrix);
        //}
        //else if (4 == num)
        //{
        //    int number = 1;
        //    cout << "���������" << endl;
        //    Matrix m = CreateMatrix();
        //    cout << "��������ֵ" << endl;
        //    cin >> number;
        //    cout << "����Ϊ��" << endl;
        //    printMatrix(m);
        //    cout << "��ֵΪ��" << endl;
        //    cout << number << endl;
        //    printMatrix(numMul(m, number));
        //}
        //else 

        cout << "���س�����....";
        /*cout << m1.matrix;*/
        /*free(m1.matrix);
        free(m2);*/
        getchar();
        getchar();
        system("cls");
    } while (0 == num || 1 == num );
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
