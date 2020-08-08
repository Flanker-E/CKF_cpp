// CKFilter.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include <vector>        //提供向量头文件
#include <algorithm>     // 算法头文件，提供迭代器
#include <fstream>       //提供文件头文件
#include <iomanip>       //C++输出精度控制需要
#include "Matrix_Operation.h"
using namespace std;

#define GM 3.986005e+014
#define J2 1.082628e-03
#define Re 6378137
#define We 7.292115e-5
#define Copies 1

void InOrbitPredict(Matrix PositionOut, Matrix AccellerationOut, Matrix PositionIn);

int main()
{
	////////////read txt file///////////////
	std::vector<double> data;
	vector<double>::iterator it;
	ifstream datatxt;
	datatxt.open("GPS_test.txt");
	double d;
	while (datatxt >> d)
		data.push_back(d);//将数据压入堆栈。
	datatxt.close();

	////////////init data///////////////
	int sizeofdata = (int)data.size(); //size of data
	int UnitSize = 142;   // Each second we have 142 variables traced now.
	int lengthData = floor(sizeofdata / UnitSize);//how much epochs of data

	Matrix UTC = ZerosMatrix(lengthData,1); //in seconds
	Matrix ARMOutputUserPosX = ZerosMatrix(lengthData, 1); // user position ARMOutputUserPosX
	Matrix ARMOutputUserPosY = ZerosMatrix(lengthData, 1);
	Matrix ARMOutputUserPosZ = ZerosMatrix(lengthData, 1);
	Matrix ARMOutputUserPosDt = ZerosMatrix(lengthData, 1);
	Matrix ARMOutputUserVelocityX = ZerosMatrix(lengthData, 1); // user velocity
	Matrix ARMOutputUserVelocityY = ZerosMatrix(lengthData, 1);
	Matrix ARMOutputUserVelocityZ = ZerosMatrix(lengthData, 1);
	Matrix State = ZerosMatrix(lengthData, 1);

	Matrix PseudoRange = ZerosMatrix(lengthData, 12);
	Matrix NumberGPSSpaceVehiclesInView = ZerosMatrix(lengthData, 1);
	Matrix GPSSpaceVehiclePosX = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclePosZ = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclePosY = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesVelocityX = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesVelocityZ = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesVelocityY = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesdL1 = ZerosMatrix(lengthData, 12);

	///////////process data///////////
	for (int i = 0; i < lengthData; i++)
	{
		UTC.matrix[i][0] = *(data.begin() + i * UnitSize + 1 - 1);
		ARMOutputUserPosX.matrix[i][0] = *(data.begin() + i * UnitSize + 2 - 1);
		ARMOutputUserPosY.matrix[i][0] = *(data.begin() + i * UnitSize + 3 - 1);
		ARMOutputUserPosZ.matrix[i][0] = *(data.begin() + i * UnitSize + 4 - 1);
		ARMOutputUserVelocityX.matrix[i][0] = *(data.begin() + i * UnitSize + 5 - 1) / 256;
		ARMOutputUserVelocityY.matrix[i][0] = *(data.begin() + i * UnitSize + 6 - 1) / 256;
		ARMOutputUserVelocityZ.matrix[i][0] = *(data.begin() + i * UnitSize + 7 - 1) / 256;
		//cout << ARMOutputUserPosX.matrix[i][0] << endl<< ARMOutputUserPosY.matrix[i][0]<<endl;
		State.matrix[i][0] = *(data.begin() + i * UnitSize + 9 - 1);
		int elementNumber = 3;
		int UsedSvNumloop = ((int)State.matrix[i][0]- (int)State.matrix[i][0] % 16)/16;
		for (int ii = 0; ii < UsedSvNumloop; ii++)
			if (*(data.begin() + i * UnitSize + 9 + ii * elementNumber + 2 - 1) > 0)
			{
				//cout << *(data.begin() + i * UnitSize + 9 + ii * elementNumber + 2 - 1) << endl;
				NumberGPSSpaceVehiclesInView.matrix[i][0]++;
				PseudoRange.matrix[i][(int)NumberGPSSpaceVehiclesInView.matrix[i][0]] = *(data.begin() + i * UnitSize + 9 + ii * elementNumber + 3 - 1);
			}
		//cout << UsedSvNumloop<<endl<<State.matrix[i][0] << endl<< NumberGPSSpaceVehiclesInView.matrix[i][0]<<endl;
		int elementNumber2 = 8;
		for (int j = 0; j < NumberGPSSpaceVehiclesInView.matrix[i][0]; j++)
		{
			GPSSpaceVehiclePosX.matrix[i][j] = *(data.begin() + i * UnitSize + j * elementNumber2 + 1 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclePosY.matrix[i][j] = *(data.begin() + i * UnitSize + j * elementNumber2 + 2 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclePosZ.matrix[i][j] = *(data.begin() + i * UnitSize + j * elementNumber2 + 3 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclesdL1.matrix[i][j] = *(data.begin() + i * UnitSize + j * elementNumber2 + 4 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclesVelocityX.matrix[i][j] = *(data.begin() + i * UnitSize + j * elementNumber2 + 5 + 10 + elementNumber * 12 -1) / 256;
			GPSSpaceVehiclesVelocityY.matrix[i][j] = *(data.begin() + i * UnitSize + j * elementNumber2 + 6 + 10 + elementNumber * 12 -1) / 256;
			GPSSpaceVehiclesVelocityZ.matrix[i][j] = *(data.begin() + i * UnitSize + j * elementNumber2 + 7 + 10 + elementNumber * 12 -1) / 256;
			//cout << GPSSpaceVehiclesdL1.matrix[i][j] << endl;
		}
	}

	////////////////adaptive filter////////////////
	Matrix PseudoRange2 = ZerosMatrix(lengthData, 12);
	Matrix NumberGPSSpaceVehiclesInView2 = ZerosMatrix(lengthData, 1);
	Matrix GPSSpaceVehiclePosX2 = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclePosZ2 = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclePosY2 = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesVelocityX2 = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesVelocityZ2 = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesVelocityY2 = ZerosMatrix(lengthData, 12);
	Matrix GPSSpaceVehiclesdL12 = ZerosMatrix(lengthData, 12);

	Matrix lmsXu = ZerosMatrix(6, lengthData); //x, y, z, vx, vy, vz
	Matrix XuKKfinal = ZerosMatrix(6, lengthData);
	Matrix XuKKtime = ZerosMatrix(6, lengthData);
	Matrix XuPredicted = ZerosMatrix(6, lengthData);
	Matrix XuKKlaststep = ZerosMatrix(6, 1);
	Matrix lmsXuErrorPower = ZerosMatrix(6, lengthData);
	Matrix P_K_Kfinal = ZerosMatrix(6, lengthData);
	Matrix P_K_Ktime = ZerosMatrix(6, lengthData);
	Matrix P_K_Klaststep = ZerosMatrix(6, lengthData);
	Matrix PosFinal = ZerosMatrix(6, lengthData);
	float Ppower = 0.02*0.02;
	float Vpower = 0.01*0.01;
	Matrix DelatYRangeAll = ZerosMatrix(6, lengthData);
	Matrix AccellerationOut = ZerosMatrix(3, lengthData);
	Matrix AccellerationOutNoUse = ZerosMatrix(3, 1);
	Matrix Q = ZerosMatrix(6, 6); //obzervation Matrix
	Q.matrix[0][0] = Ppower; Q.matrix[1][1] = Ppower; Q.matrix[2][2] = Ppower; Q.matrix[3][3] = Vpower; Q.matrix[4][4] = Vpower; Q.matrix[5][5] = Ppower;
	Matrix R = ZerosMatrix(6, 6);
	passMatrix(R, Q);
	passMatrix(P_K_Kfinal, Q);
	//R.matrix[0][0] = Ppower; R.matrix[1][1] = Ppower; R.matrix[2][2] = Ppower; R.matrix[3][3] = Vpower; R.matrix[4][4] = Vpower; R.matrix[5][5] = Ppower;
	//P_K_Kfinal = [Ppower, 0, 0, 0, 0, 0; 0, Ppower, 0, 0, 0, 0; 0, 0, Ppower, 0, 0, 0; 0, 0, 0, Vpower, 0, 0; 0, 0, 0, 0, Vpower, 0; 0, 0, 0, 0, 0, Ppower]; %obzervation Matrix




	////////////////print txt file///////////
	//int i = 0;
	ofstream out("out.txt");
	for (it = data.begin(); it != data.end(); it++)
	{
		out << setprecision(16) << *it << endl;
		//cout << "V[" << i << "]=" << setprecision(16) << *it << endl;
		//i++;
	}
	out.close();
	return 0;
}

void InOrbitPredict(Matrix PositionOut, Matrix AccellerationOut, Matrix PositionIn)
{
	Matrix New_Position = ZerosMatrix(9, 1);
	passMatrix(New_Position, PositionIn);
	double r = sqrt(pow(New_Position.matrix[0][0], 2) + pow(New_Position.matrix[1][0], 2) + pow(New_Position.matrix[2][0], 2));
	New_Position.matrix[6][0] = -(GM / pow(r, 3)) * New_Position.matrix[0][0] * (1 + 1.5 * J2 * pow((Re / r), 2) * (1 - 5 * pow((New_Position.matrix[2][0] / r), 2))) + pow(We, 2) * New_Position.matrix[0][0] + 2 * We * New_Position.matrix[4][0];
	New_Position.matrix[7][0] = -(GM / pow(r, 3)) * New_Position.matrix[1][0] * (1 + 1.5 * J2 * pow((Re / r), 2) * (1 - 5 * pow((New_Position.matrix[2][0] / r), 2))) + pow(We, 2) * New_Position.matrix[1][0] - 2 * We * New_Position.matrix[3][0];
	New_Position.matrix[8][0] = -(GM / pow(r, 3)) * New_Position.matrix[2][0] * (1 + 1.5 * J2 * pow((Re / r), 2) * (3 - 5 * pow((New_Position.matrix[2][0] / r), 2)));
	for (int i = 0; i < 3; i++)
	{
		New_Position.matrix[i][0] = New_Position.matrix[i][0] + New_Position.matrix[i + 3][0] / pow(Copies, 2) + 0.5*New_Position.matrix[i + 6][0] / pow(Copies, 2); //位置
		New_Position.matrix[i + 3][0] = New_Position.matrix[i + 3][0] + New_Position.matrix[i + 6][0] / pow(Copies, 2); //速度
	}
	passMatrix(PositionOut, New_Position, 1, 0, 1, 0, 1, 6);
	passMatrix(AccellerationOut, New_Position, 1, 0, 1, 0, 7, 9);
}
