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
	datatxt.open("GPSTrace.txt");
	double d;
	while (datatxt >> d)
		data.push_back(d);//将数据压入堆栈。
	datatxt.close();

	std::vector<double> TruePosInfodata;
	//vector<double>::iterator it;
	//ifstream datatxt;
	datatxt.open("TrueData.TXT");
	//double d;
	while (datatxt >> d)
		TruePosInfodata.push_back(d);//将数据压入堆栈。
	datatxt.close();

	////////////init data///////////////
	int sizeofdata = (int)data.size(); //size of data
	int UnitSize = 142;   // Each second we have 142 variables traced now.
	int lengthData = floor(sizeofdata / UnitSize);//how much epochs of data

	int sizeofTruePosInfodata = (int)TruePosInfodata.size();
	int lengthTruePosInfodata = sizeofTruePosInfodata / 20;
	Matrix TruePosInfo = ZerosMatrix(6, lengthData);//6, lengthData
	it = TruePosInfodata.begin();
	for (int i = 0; i < lengthData; i++)
	{
		it += 2;
		for (int ii = 0; ii < 6; ii++)
		{
			TruePosInfo.matrix[ii][i] = *it;
			it++;
		}
		it += 12;
	}
	Matrix UTC = ZerosMatrix(1,lengthData); //1,lengthData in seconds
	Matrix ARMOutputUserPosX = ZerosMatrix(1,lengthData); //1,lengthData user position ARMOutputUserPosX
	Matrix ARMOutputUserPosY = ZerosMatrix(1,lengthData);//1,lengthData
	Matrix ARMOutputUserPosZ = ZerosMatrix(1,lengthData);//1,lengthData
	Matrix ARMOutputUserPosDt = ZerosMatrix(1,lengthData);//1,lengthData
	Matrix ARMOutputUserVelocityX = ZerosMatrix(1,lengthData); //1,lengthData user velocity
	Matrix ARMOutputUserVelocityY = ZerosMatrix(1,lengthData);//1,lengthData
	Matrix ARMOutputUserVelocityZ = ZerosMatrix(1,lengthData);//1,lengthData
	Matrix State = ZerosMatrix(1,lengthData);//1,lengthData

	Matrix PseudoRange = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix NumberGPSSpaceVehiclesInView = ZerosMatrix(1,lengthData);//1,lengthData
	Matrix GPSSpaceVehiclePosX = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclePosZ = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclePosY = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesVelocityX = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesVelocityZ = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesVelocityY = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesdL1 = ZerosMatrix(12, lengthData);//12,lengthData

	///////////process data///////////
	for (int i = 0; i < lengthData; i++)
	{
		UTC.matrix[0][i] = *(data.begin() + i * UnitSize + 1 - 1);
		ARMOutputUserPosX.matrix[0][i] = *(data.begin() + i * UnitSize + 2 - 1);
		ARMOutputUserPosY.matrix[0][i] = *(data.begin() + i * UnitSize + 3 - 1);
		ARMOutputUserPosZ.matrix[0][i] = *(data.begin() + i * UnitSize + 4 - 1);
		ARMOutputUserVelocityX.matrix[0][i] = *(data.begin() + i * UnitSize + 5 - 1) / 256;
		ARMOutputUserVelocityY.matrix[0][i] = *(data.begin() + i * UnitSize + 6 - 1) / 256;
		ARMOutputUserVelocityZ.matrix[0][i] = *(data.begin() + i * UnitSize + 7 - 1) / 256;
		//cout << ARMOutputUserPosX.matrix[0][i] << endl<< ARMOutputUserPosY.matrix[0][i]<<endl;
		State.matrix[0][i] = *(data.begin() + i * UnitSize + 9 - 1);
		int elementNumber = 3;
		int UsedSvNumloop = ((int)State.matrix[0][i]- (int)State.matrix[0][i] % 16)/16;
		for (int ii = 0; ii < UsedSvNumloop; ii++)
			if (*(data.begin() + i * UnitSize + 9 + ii * elementNumber + 2 - 1) > 0)
			{
				//cout << *(data.begin() + i * UnitSize + 9 + ii * elementNumber + 2 - 1) << endl;
				NumberGPSSpaceVehiclesInView.matrix[0][i]++;
				PseudoRange.matrix[(int)NumberGPSSpaceVehiclesInView.matrix[0][i]-1][i] = *(data.begin() + i * UnitSize + 9 + ii * elementNumber + 3 - 1);
			}
		//cout << UsedSvNumloop<<endl<<State.matrix[0][i] << endl<< NumberGPSSpaceVehiclesInView.matrix[0][i]<<endl;
		int elementNumber2 = 8;
		for (int j = 0; j < NumberGPSSpaceVehiclesInView.matrix[0][i]; j++)
		{
			GPSSpaceVehiclePosX.matrix[j][i] = *(data.begin() + i * UnitSize + j * elementNumber2 + 1 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclePosY.matrix[j][i] = *(data.begin() + i * UnitSize + j * elementNumber2 + 2 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclePosZ.matrix[j][i] = *(data.begin() + i * UnitSize + j * elementNumber2 + 3 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclesdL1.matrix[j][i] = *(data.begin() + i * UnitSize + j * elementNumber2 + 4 + 10 + elementNumber * 12 -1);
			GPSSpaceVehiclesVelocityX.matrix[j][i] = *(data.begin() + i * UnitSize + j * elementNumber2 + 5 + 10 + elementNumber * 12 -1) / 256;
			GPSSpaceVehiclesVelocityY.matrix[j][i] = *(data.begin() + i * UnitSize + j * elementNumber2 + 6 + 10 + elementNumber * 12 -1) / 256;
			GPSSpaceVehiclesVelocityZ.matrix[j][i] = *(data.begin() + i * UnitSize + j * elementNumber2 + 7 + 10 + elementNumber * 12 -1) / 256;
			//cout << GPSSpaceVehiclesdL1.matrix[j][i] << endl;
		}
	}

	////////////////adaptive filter////////////////
	Matrix PseudoRange2 = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix NumberGPSSpaceVehiclesInView2 = ZerosMatrix(1,lengthData);//1,lengthData
	Matrix GPSSpaceVehiclePosX2 = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclePosZ2 = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclePosY2 = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesVelocityX2 = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesVelocityZ2 = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesVelocityY2 = ZerosMatrix(12, lengthData);//12,lengthData
	Matrix GPSSpaceVehiclesdL12 = ZerosMatrix(12, lengthData);//12,lengthData

	Matrix lmsXu = ZerosMatrix(6, lengthData); //6,lengthData,x, y, z, vx, vy, vz
	Matrix XuKKfinal = ZerosMatrix(6, lengthData);//6,lengthData
	Matrix XuKKtime = ZerosMatrix(6, 1);//6,1
	Matrix XuPredicted = ZerosMatrix(6, lengthData);//6,lengthData
	Matrix XuKKlaststep = ZerosMatrix(6, 1);//6,1
	Matrix lmsXuErrorPower = ZerosMatrix(6, lengthData);//6,lengthData
	Matrix P_K_Kfinal = ZerosMatrix(6, 6);//6,6 symmetric matrix
	Matrix P_K_Ktime = ZerosMatrix(6, 6);//6,6
	Matrix P_K_Klaststep = ZerosMatrix(6, 6);//6,6
	Matrix PosFinal = ZerosMatrix(6, lengthData);//6,lengthData
	double Ppower = 0.02*0.02; 
	double Vpower = 0.01*0.01;
	Matrix DeltaYRangeAll = ZerosMatrix(6, lengthData);//6,lengthData
	Matrix AccellerationOut = ZerosMatrix(3, lengthData);//3,lengthData
	Matrix Q = ZerosMatrix(6, 6); //6,6,obzervation Matrix
	Q.matrix[0][0] = Ppower; Q.matrix[1][1] = Ppower; Q.matrix[2][2] = Ppower; Q.matrix[3][3] = Vpower; Q.matrix[4][4] = Vpower; Q.matrix[5][5] = Ppower;
	Ppower = 1*1;
	Vpower = 0.02*0.02;
	Matrix R = ZerosMatrix(6, 6);//6,6
	R.matrix[0][0] = Ppower; R.matrix[1][1] = Ppower; R.matrix[2][2] = Ppower; R.matrix[3][3] = Vpower; R.matrix[4][4] = Vpower; R.matrix[5][5] = Ppower;
	Matrix PredictedPosInfo = ZerosMatrix(6, 1);//6,1
	passMatrix(P_K_Kfinal, R);
	//printMatrix(R);	printMatrix(P_K_Kfinal);
	for (int i = 0; i < lengthData; i++)
	{
		lmsXu.matrix[0][i] = ARMOutputUserPosX.matrix[0][i];
		lmsXu.matrix[1][i] = ARMOutputUserPosY.matrix[0][i];
		lmsXu.matrix[2][i] = ARMOutputUserPosZ.matrix[0][i];
		lmsXu.matrix[3][i] = ARMOutputUserVelocityX.matrix[0][i];
		lmsXu.matrix[4][i] = ARMOutputUserVelocityY.matrix[0][i];
		lmsXu.matrix[5][i] = ARMOutputUserVelocityZ.matrix[0][i];
	}
	//printMatrix(lmsXu);
	passMatrix(XuKKfinal, lmsXu);	
	passMatrix(XuPredicted, lmsXu);
	//why start at 100?
	int LoopStart = 0;// 100;
	int SvNum = NumberGPSSpaceVehiclesInView.matrix[0][LoopStart];
	Matrix Doppler = ZerosMatrix(12, 1);//12,1
	Matrix DopplerError = ZerosMatrix(12, 1);//12,1
	for (int i = 0; i < SvNum; i++)
	{
		double dx = ARMOutputUserPosX.matrix[0][LoopStart] - GPSSpaceVehiclePosX.matrix[i][LoopStart];
		double dy = ARMOutputUserPosY.matrix[0][LoopStart] - GPSSpaceVehiclePosY.matrix[i][LoopStart];
		double dz = ARMOutputUserPosZ.matrix[0][LoopStart] - GPSSpaceVehiclePosZ.matrix[i][LoopStart];
		double dvx = ARMOutputUserVelocityX.matrix[0][LoopStart] - GPSSpaceVehiclesVelocityX.matrix[i][LoopStart];
		double dvy = ARMOutputUserVelocityY.matrix[0][LoopStart] - GPSSpaceVehiclesVelocityY.matrix[i][LoopStart];
		double dvz = ARMOutputUserVelocityZ.matrix[0][LoopStart] - GPSSpaceVehiclesVelocityZ.matrix[i][LoopStart];
		double r = sqrt(dx*dx + dy * dy + dz * dz);
		Doppler.matrix[i][0] = (dx*dvx + dy * dvy + dz * dvz) / r;
		DopplerError.matrix[i][0] = GPSSpaceVehiclesdL1.matrix[i][LoopStart] - Doppler.matrix[i][0];
	}
	double Df = meanMatrix(DopplerError, 1, SvNum);
	passMatrix(XuKKlaststep, lmsXu, 1, 0, 1, 0, 1, 0, 1, 1);
	//cout << Df << endl;
	//printMatrix(XuKKlaststep);

	///////////cubaTrue时间更新初始化////////////
	//Initialize Epsilon;
	int StateNumber = 6;//set n=6
	int m = StateNumber * 2;//m=2n=12
	Matrix CPRcv = ZerosMatrix(StateNumber, m);//6,12 StateNumber, m CPRcv=X^*_{i,k|k-1}
	//DopplerError(1:SvNum);
	Matrix Epsilon = ZerosMatrix(6, 12);//6,12
	Matrix TempEye = EyeMatrix(6);//6,6
	passMatrix(Epsilon, TempEye);
	passMatrix(Epsilon, numMul(TempEye, -1, 1), 1, 0, 6 + 1);
	Epsilon = numMul(Epsilon, (double)sqrt(6.0),1);
	//printMatrix(Epsilon);
	//////////declaration in main loop////////////
	Matrix TempA, TempB, TempC, TempOne, TempZero;//temporary matrix
	Matrix TRange, TDoppler;//SvNum,1
	Matrix TempPseu, TempTestDtMedian, TestDeltaYRange;
	Matrix Tempsdl1, TempTestDfMedian, TestDeltaYDoppler;
	double TestDtMedian, TestDfMedian;
	Matrix D;//SvNum - 1, SvNum
	Matrix D1;//2 * (SvNum - 1), 2 * SvNum
	Matrix R1;//2 * (SvNum - 1), 2 * (SvNum - 1)
	Matrix Z;//2 * (SvNum - 1), 1
	Matrix S_K_Ktime, S_K_Kmeas;//6,6
	Matrix CP, CPmeas;//6,12
	Matrix AccellerationOutNoUse;//3,1
	Matrix CPYRange, CPYDoppler;//SvNum,1
	Matrix YRange, YDoppler, DeltaYRange;//SvNum, 1
	Matrix DeltaYDoppler;
	Matrix YRobust, Zkmeas;//2 * SvNum, 1
	double DtMedian, DfMedian, RangeThreshSmall, DopplerThreshSmall;
	Matrix CPZmeas, CPYmeas, CPY;//2* SvNum, m
	Matrix Pzz_K_Kmeas, Pxz_K_Kmeas;
	Matrix Wk;
	double alpha;
	Matrix FinalRange, FinalDoppler;//SvNum,1
	double DtFinal, DopplerFinal, RangePower, PowerYRange, PowerYDopper, ProporgatedPowerYRange, ProporgatedPowerYDopper;
	Matrix R21, R22;//2 * SvNum, 2 * SvNum
	Matrix diagR21, diagR22;//2 * SvNum, 1

	///////////////////////////main loop//////////////////////////////////
	for (int i = LoopStart + 1; i < lengthData; i++)
	{
		/*if ((i % 100) == 0)*/
		//cout << i << endl;
		passMatrix(P_K_Klaststep, P_K_Kfinal);
		passMatrix(XuKKlaststep, XuKKfinal, 1, 0, 1, 0, 1, 0, i, i);
		SvNum = NumberGPSSpaceVehiclesInView.matrix[0][LoopStart];
		TempA = ZerosMatrix(6, 1);
		TempB = ZerosMatrix(3, 1);
		TempC = ZerosMatrix(6, 1);
		InOrbitPredict(PredictedPosInfo, TempB, XuKKlaststep);
		passMatrix(AccellerationOut, TempB, 1, 0, i + 1, i + 1);
		passMatrix(TempC, TruePosInfo, 1, 0, 1, 0, 1, 0, i, i);
		InOrbitPredict(TempA, TempB, TempC);
		passMatrix(XuPredicted, TempA, 1, 0, i + 1, i + 1);
		freeMatrix(TempA); freeMatrix(TempB); freeMatrix(TempC);
		//printMatrix(P_K_Klaststep);
		//printMatrix(XuKKlaststep);  printMatrix(AccellerationOut); printMatrix(XuPredicted);
		//cout << SvNum;
		TRange = ZerosMatrix(SvNum, 1);//SvNum,1
		TDoppler = ZerosMatrix(SvNum, 1);//SvNum,1
		//Matrix TDopplerLMS = ZerosMatrix(SvNum, 1);
		for (int ii = 0; ii < SvNum; ii++)
		{
			double dx = PredictedPosInfo.matrix[0][0] - GPSSpaceVehiclePosX.matrix[ii][i];
			double dy = PredictedPosInfo.matrix[1][0] - GPSSpaceVehiclePosY.matrix[ii][i];
			double dz = PredictedPosInfo.matrix[2][0] - GPSSpaceVehiclePosZ.matrix[ii][i];
			double dvx = PredictedPosInfo.matrix[3][0] - GPSSpaceVehiclesVelocityX.matrix[ii][i];
			double dvy = PredictedPosInfo.matrix[4][0] - GPSSpaceVehiclesVelocityY.matrix[ii][i];
			double dvz = PredictedPosInfo.matrix[5][0] - GPSSpaceVehiclesVelocityZ.matrix[ii][i];
			double r = sqrt(dx* dx + dy * dy + dz * dz);
			TRange.matrix[ii][0] = r;
			TDoppler.matrix[ii][0] = (dx*dvx + dy * dvy + dz * dvz) / r;
			//DopplerError.matrix[i][0] = GPSSpaceVehiclesdL1.matrix[i][LoopStart] - Doppler.matrix[i][0];
		}
		//printMatrix(TRange); printMatrix(TDoppler);

		TempPseu = ZerosMatrix(SvNum, 1);
		passMatrix(TempPseu, PseudoRange, 1, 0, 1, 0, 1, SvNum, i + 1, i + 1);
		//printMatrix(sub(TempPseu, TRange));
		TestDtMedian = medianMatrix(sub(TempPseu, TRange), 0, 1);
		TempTestDtMedian = OnesMatrix(SvNum, 1, TestDtMedian);
		TestDeltaYRange = sub(sub(TempPseu, TRange, 11), TempTestDtMedian, 11);
		
		Tempsdl1 = ZerosMatrix(SvNum, 1);
		passMatrix(Tempsdl1, GPSSpaceVehiclesdL1, 1, 0, 1, 0, 1, SvNum, i + 1, i + 1);
		TestDfMedian = medianMatrix(sub(Tempsdl1, TDoppler), 0, 1);
		TempTestDfMedian = OnesMatrix(SvNum, 1, TestDfMedian);
		TestDeltaYDoppler = sub(sub(Tempsdl1, TDoppler, 11), TempTestDfMedian, 11);
		//double Df = meanMatrix(DopplerError, 1, SvNum);
		

		int index = 0;
		int SvToDelete = 0;
		for (int ii = 0; ii < SvNum; ii++)
		{
			if (abs(TestDeltaYRange.matrix[ii][0]) < 10 && abs(TestDeltaYDoppler.matrix[ii][0]) < 0.2)
			{
				GPSSpaceVehiclePosX2.matrix[index][i] = GPSSpaceVehiclePosX.matrix[ii][i];
				GPSSpaceVehiclePosY2.matrix[index][i] = GPSSpaceVehiclePosY.matrix[ii][i];
				GPSSpaceVehiclePosZ2.matrix[index][i] = GPSSpaceVehiclePosZ.matrix[ii][i];
				PseudoRange2.matrix[index][i] = PseudoRange.matrix[ii][i];
				GPSSpaceVehiclesVelocityX2.matrix[index][i] = GPSSpaceVehiclesVelocityX.matrix[ii][i];
				GPSSpaceVehiclesVelocityY2.matrix[index][i] = GPSSpaceVehiclesVelocityY.matrix[ii][i];
				GPSSpaceVehiclesVelocityZ2.matrix[index][i] = GPSSpaceVehiclesVelocityZ.matrix[ii][i];
				GPSSpaceVehiclesdL12.matrix[index][i] = GPSSpaceVehiclesdL1.matrix[ii][i];
				index++;
			}
			else
				SvToDelete++;
		}
		freeMatrix(TestDeltaYRange); freeMatrix(TestDeltaYDoppler);
		//update available device number
		SvNum -= SvToDelete;
		//cout << SvNum << endl;
		//printMatrix(GPSSpaceVehiclePosX2); printMatrix(PseudoRange2); printMatrix(GPSSpaceVehiclesdL12);

		D = ZerosMatrix(SvNum - 1, SvNum);
		TempEye = EyeMatrix(SvNum - 1);
		TempOne = OnesMatrix(SvNum - 1, 1, -1);
		passMatrix(D, TempEye);
		passMatrix(D, TempOne, 1, 0, SvNum, SvNum);
		freeMatrix(TempEye);
		freeMatrix(TempOne);
		TempZero = ZerosMatrix(SvNum - 1, SvNum);
		D1 = ZerosMatrix(2 * (SvNum - 1), 2 * SvNum);
		passMatrix(D1, D);
		passMatrix(D1, TempZero, 1, 0, SvNum + 1);
		passMatrix(D1, TempZero, SvNum, 0, 1, 0);
		passMatrix(D1, D, SvNum, 0, SvNum + 1);
		freeMatrix(D);
		freeMatrix(TempZero);
		TempOne = OnesMatrix(1, 2 * SvNum, Vpower);
		for (int ii = 0; ii < SvNum; ii++)
			TempOne.matrix[0][ii] = Ppower;
		R = diagMatrix(TempOne,1);
		R1 = Mul(Mul(D1, R, 01), TransposeMatrix(D1), 11);
		TempOne = ZerosMatrix(2*SvNum, 1);
		passMatrix(TempOne, PseudoRange2, 1, 0, 1, 0, 1, SvNum, i + 1, i + 1);
		//printMatrix(TempOne); printMatrix(D1); printMatrix(GPSSpaceVehiclesdL12);
		passMatrix(TempOne, GPSSpaceVehiclesdL12, SvNum + 1, 0, 1, 0, 1, SvNum, i + 1, i + 1);
		Z = Mul(D1, TempOne, 01);
		freeMatrix(Z);
		//freeMatrix(TempOne);
		//printMatrix(Z);

		////////////////cubaTrue时间更新 Next we goes to time update.

		S_K_Ktime = Cholesky(P_K_Klaststep);
		CP = add(Mul(S_K_Ktime, Epsilon,10), repmatMatrix(XuKKlaststep, 1, 12), 11); //6,12
		//将上一时刻CP点递推到下一时刻，传播cubaTrue点
		//CPRcv = ZerosMatrix(6, 2 * StateNumber);//6, 2 * StateNumber
		for (int j = 0; j < 2 * StateNumber; j++)
		{
			TempOne = ZerosMatrix(6, 1);
			AccellerationOutNoUse = ZerosMatrix(3, 1);
			passMatrix(TempOne, CP, 1, 0, 1, 0, 1, 0, j + 1, j + 1);
			InOrbitPredict(TempOne, AccellerationOutNoUse, TempOne);
			passMatrix(CPRcv, TempOne, 1, 0, j + 1, j + 1);
			freeMatrix(TempOne);
			freeMatrix(AccellerationOutNoUse);
		}
		freeMatrix(CP);
		for (int ii = 0; ii < 6; ii++)
			XuKKtime.matrix[ii][0] = meanMatrix(CPRcv, ii + 1, ii + 1);
		//TempA = CPRcv;
		P_K_Ktime = add(Q,numMul(Mul(sub(CPRcv, repmatMatrix(XuKKtime, 1, 12),01),TransposeMatrix(sub(CPRcv, repmatMatrix(XuKKtime, 1, 12),01),1),11), 1.0 / m,1) ,01);//估计 k 时刻的状态误差协方差预测值
		//printMatrix(P_K_Ktime);

		///////////测量更新 meas update///////////

		S_K_Kmeas = Cholesky(P_K_Ktime);//6,6
		CPmeas = add(Mul(S_K_Kmeas, Epsilon,10), repmatMatrix(XuKKtime, 1, 12), 11); //6,12	
		//M estimator to constraint outliers.
		CPYRange = ZerosMatrix(SvNum, m);//SvNum,1
		CPYDoppler = ZerosMatrix(SvNum, m);//SvNum,1
		for (int ii = 0; ii < SvNum; ii++)
		{
			for (int j = 0; j < 2 * StateNumber; j++)
			{
				double dx = CPmeas.matrix[0][j] - GPSSpaceVehiclePosX2.matrix[ii][i];
				double dy = CPmeas.matrix[1][j] - GPSSpaceVehiclePosY2.matrix[ii][i];
				double dz = CPmeas.matrix[2][j] - GPSSpaceVehiclePosZ2.matrix[ii][i];
				double dvx = CPmeas.matrix[3][j] - GPSSpaceVehiclesVelocityX2.matrix[ii][i];
				double dvy = CPmeas.matrix[4][j] - GPSSpaceVehiclesVelocityY2.matrix[ii][i];
				double dvz = CPmeas.matrix[5][j] - GPSSpaceVehiclesVelocityZ2.matrix[ii][i];
				double r = sqrt(dx* dx + dy * dy + dz * dz);
				CPYRange.matrix[ii][j] = r;
				CPYDoppler.matrix[ii][j] = (dx*dvx + dy * dvy + dz * dvz) / r;
				//DopplerError.matrix[i][0] = GPSSpaceVehiclesdL1.matrix[i][LoopStart] - Doppler.matrix[i][0];
			}
		}
		YRange = ZerosMatrix(SvNum, 1);
		YDoppler = ZerosMatrix(SvNum, 1);
		for (int ii = 0; ii < SvNum; ii++)
			YRange.matrix[ii][0] = meanMatrix(CPYRange, ii + 1, ii + 1);
		for (int ii = 0; ii < SvNum; ii++)
			YDoppler.matrix[ii][0] = meanMatrix(CPYDoppler, ii + 1, ii + 1);
		//printMatrix(TRange); printMatrix(TDoppler);
		TempPseu = ZerosMatrix(SvNum, 1);
		passMatrix(TempPseu, PseudoRange2, 1, 0, 1, 0, 1, SvNum, i + 1, i + 1);
		DtMedian = medianMatrix(sub(TempPseu, YRange), 0, 1);
		DeltaYRange = sub(sub(TempPseu, YRange, 10), OnesMatrix(SvNum, 1, DtMedian), 11);

		Tempsdl1 = ZerosMatrix(SvNum, 1);
		passMatrix(Tempsdl1, GPSSpaceVehiclesdL12, 1, 0, 1, 0, 1, SvNum, i + 1, i + 1);
		DfMedian = medianMatrix(sub(Tempsdl1, YDoppler), 0, 1);
		DeltaYDoppler = sub(sub(Tempsdl1, YDoppler, 10), OnesMatrix(SvNum, 1, DfMedian), 11);
		//double Df = meanMatrix(DopplerError, 1, SvNum);
		YRobust = ZerosMatrix(2 * SvNum, 1);
		RangeThreshSmall = 0.8; //sqrt(Ppower) * 3; 待查证？
		DopplerThreshSmall = sqrt(Vpower) * 3;

		////////////////影响函数////////////////
		for (int ii = 0; ii < SvNum; ii++)
		{
			if(DeltaYRange.matrix[ii][0]>RangeThreshSmall)
				YRobust.matrix[ii][0] = YRange.matrix[ii][0] + RangeThreshSmall;
			else if(DeltaYRange.matrix[ii][0] < -RangeThreshSmall)
				YRobust.matrix[ii][0] = YRange.matrix[ii][0] - RangeThreshSmall;
			else
				YRobust.matrix[ii][0] = YRange.matrix[ii][0] + DeltaYRange.matrix[ii][0];

			if (DeltaYDoppler.matrix[ii][0] > DopplerThreshSmall)
				YRobust.matrix[ii+SvNum][0] = YDoppler.matrix[ii][0] + DopplerThreshSmall;
			else if(DeltaYDoppler.matrix[ii][0] < -DopplerThreshSmall)
				YRobust.matrix[ii + SvNum][0] = YDoppler.matrix[ii][0] - DopplerThreshSmall;
			else
				YRobust.matrix[ii + SvNum][0] = YDoppler.matrix[ii][0] + DeltaYDoppler.matrix[ii][0];
		}

		//估计相关协方差阵
		TempOne = ZerosMatrix(2*SvNum, m);
		passMatrix(TempOne, CPYRange, 1, 0, 1, 0, 1, SvNum);
		freeMatrix(CPYRange);
		//printMatrix(TempOne); printMatrix(D1); printMatrix(GPSSpaceVehiclesdL12);
		passMatrix(TempOne, CPYDoppler, SvNum + 1, 0, 1, 0, 1, SvNum);
		freeMatrix(CPYDoppler);

		CPYmeas = TempOne;
		CPZmeas = Mul(D1, TempOne);
		
		TempOne = ZerosMatrix(2 * SvNum, 1);
		passMatrix(TempOne, YRange, 1, 0, 1, 0, 1, SvNum);
		passMatrix(TempOne, YDoppler, SvNum + 1, 0, 1, 0, 1, SvNum);
		freeMatrix(YRange); freeMatrix(YDoppler);
		CPY = TempOne;
		Zkmeas = Mul(D1, TempOne);

		Pzz_K_Kmeas = add(R1, numMul(Mul(sub(CPZmeas, repmatMatrix(Zkmeas, 1, m), 01), TransposeMatrix(sub(CPZmeas, repmatMatrix(Zkmeas, 1, m), 01), 1), 11), 1.0 / m, 1), 11);
		Pxz_K_Kmeas = numMul(Mul(sub(CPmeas, repmatMatrix(XuKKtime, 1, m), 11), TransposeMatrix(sub(CPZmeas, repmatMatrix(Zkmeas, 1, m), 11), 1), 11), 1.0 / m, 1);
		//估计卡尔曼增益
		Wk = Mul(Pxz_K_Kmeas, CholeskyInverse(Pzz_K_Kmeas), 01);
		//K 时刻状态估计值
		passMatrix(XuKKfinal, add(XuKKtime, Mul(Wk, sub(Mul(D1, YRobust,10), Zkmeas, 10), 01), 01), 1, 0, i + 1, i + 1);
		//K 时刻状态误差协方差估计值
		passMatrix(P_K_Kfinal, sub(P_K_Ktime, Mul(Wk, Mul(Pzz_K_Kmeas, TransposeMatrix(Wk), 01), 01), 11));


		alpha = -0.0;
		TempA = ZerosMatrix(6, 1);
		passMatrix(TempA, XuKKfinal, 1, 0, 1, 0, 1, 0, i + 1, i + 1);
		TempB = ZerosMatrix(6, 1);
		passMatrix(TempB, XuKKfinal, 1, 0, 1, 0, 1, 0, i, i);
		TempC = ZerosMatrix(6, 1);
		passMatrix(TempC, TruePosInfo, 1, 0, 1, 0, 1, 0, i + 1, i + 1);
		TempOne = ZerosMatrix(6, 1);
		passMatrix(TempOne, AccellerationOut, 1, 0, 1, 0, 1, 0, i + 1, i + 1);
		//printMatrix(add(numMul(TempA, (1 - alpha)), numMul(TempB, alpha)));
		//cout << Ppower << "   " << Vpower << endl;
		//printMatrix(add(TempC, numMul(TempOne, 0.5)));
		TempOne = sub(add(numMul(TempA, (1 - alpha), 1), numMul(TempB, alpha, 1), 11), add(TempC, numMul(TempOne, 0.5, 1), 11), 11);
		passMatrix(PosFinal, TempOne, 1, 0, i + 1, i + 1);
		freeMatrix(TempOne);
		

		FinalRange = ZerosMatrix(SvNum, 1);
		FinalDoppler = ZerosMatrix(SvNum, 1);//SvNum,1
		for (int ii = 0; ii < SvNum; ii++)
		{
			double dx = XuKKfinal.matrix[0][i] - GPSSpaceVehiclePosX2.matrix[ii][i];
			double dy = XuKKfinal.matrix[1][i] - GPSSpaceVehiclePosY2.matrix[ii][i];
			double dz = XuKKfinal.matrix[2][i] - GPSSpaceVehiclePosZ2.matrix[ii][i];
			double dvx = XuKKfinal.matrix[3][i] - GPSSpaceVehiclesVelocityX2.matrix[ii][i];
			double dvy = XuKKfinal.matrix[4][i] - GPSSpaceVehiclesVelocityY2.matrix[ii][i];
			double dvz = XuKKfinal.matrix[5][i] - GPSSpaceVehiclesVelocityZ2.matrix[ii][i];
			double r = sqrt(dx* dx + dy * dy + dz * dz);
			FinalRange.matrix[ii][0] = r;
			FinalDoppler.matrix[ii][0] = (dx*dvx + dy * dvy + dz * dvz) / r;//.后12位
			//DopplerError.matrix[i][0] = GPSSpaceVehiclesdL1.matrix[i][LoopStart] - Doppler.matrix[i][0];
		}

		TempA = ZerosMatrix(SvNum, 1);
		passMatrix(TempA, YRobust,1,0,1,0,1,SvNum);
		DeltaYRange = sub(TempA, FinalRange,10);
		DtFinal = meanMatrix(DeltaYRange, 1, 0, 1, 0) + DtMedian;
		//printMatrix(DeltaYRange);
		//printMatrix(DeltaYDoppler);
		TempA = ZerosMatrix(SvNum, 1);
		passMatrix(TempA, YRobust, 1, 0, 1, 0, 1 + SvNum, 2 * SvNum);
		freeMatrix(YRobust);
		DeltaYDoppler = sub(TempA, FinalDoppler, 10);
		DopplerFinal = meanMatrix(DeltaYDoppler, 1, 0, 1, 0);
		RangePower = sqrt(meanMatrix(Mul(TransposeMatrix(DeltaYRange), DeltaYRange, 10), 1, 0, 1, 0, 1) / SvNum);

		R21 = ZerosMatrix(2 * SvNum, 2 * SvNum);
		TempA = Mul(DeltaYRange, TransposeMatrix(DeltaYRange), 01);
		passMatrix(R21, TempA); 
		freeMatrix(TempA); freeMatrix(DeltaYRange);
		TempZero = ZerosMatrix(SvNum, SvNum);
		passMatrix(R21, TempZero, 1, 0, SvNum + 1);
		passMatrix(R21, TempZero, SvNum + 1, 0, 1, 0);
		TempA = Mul(DeltaYDoppler, TransposeMatrix(DeltaYDoppler), 01);
		passMatrix(R21, TempA, SvNum + 1, 0, SvNum + 1); 
		freeMatrix(TempA); freeMatrix(DeltaYDoppler);

		diagR21 = ZerosMatrix(2 * SvNum, 1);
		for (int j = 0; j < 2 * SvNum; j++)
			diagR21.matrix[j][0] = R21.matrix[j][j];
		PowerYRange = meanMatrix(diagR21, 1, SvNum);
		PowerYDopper = meanMatrix(diagR21, SvNum + 1, 2 * SvNum);
		
		//printMatrix(CPYmeas);
		//printMatrix(CPY);
		//printMatrix(repmatMatrix(CPY, 1, m));
		TempA = sub(CPYmeas, repmatMatrix(CPY, 1, m), 01);
		R22 = numMul(Mul(sub(CPYmeas, repmatMatrix(CPY, 1, m), 01), TransposeMatrix(sub(CPYmeas, repmatMatrix(CPY, 1, m), 01), 1), 11), 1.0 / m, 1);
		freeMatrix(CPYmeas);
		diagR22 = ZerosMatrix(2 * SvNum, 1);
		for (int j = 0; j < 2 * SvNum; j++)
			diagR22.matrix[j][0] = R22.matrix[j][j];
		ProporgatedPowerYRange = meanMatrix(diagR22, 1, SvNum);
		ProporgatedPowerYDopper = meanMatrix(diagR22, SvNum + 1, 2 * SvNum);
		Ppower = 0.98*Ppower + 0.02*PowerYRange;
		Vpower = 0.98*Vpower + 0.02*PowerYDopper;
	}

	////////////////print txt file///////////
	ofstream out("out.txt");
	for (int i = 0; i<lengthData; i++)
	{
		out << setprecision(16) << XuPredicted.matrix[0][i]-TruePosInfo.matrix[0][i] <<"  "
			<< XuPredicted.matrix[1][i] - TruePosInfo.matrix[1][i] << "  " << XuPredicted.matrix[2][i] - TruePosInfo.matrix[2][i] << endl;
	}
	out.close();

	//ofstream out("out.txt");
	//for (it = data.begin(); it != data.end(); it++)
	//{
	//	out << setprecision(16) << *it << endl;
	//}
	//out.close();

	//int i = 0;
	//for (it = data.begin(); it != data.end(); it++)
	//{
	//	cout << "V[" << i << "]=" << setprecision(16) << *it << endl;
	//	i++;
	//}
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
