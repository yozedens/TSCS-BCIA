#include"SetValue.h"
#include<algorithm>

//速密关系函数
double V(double u)
{
	return vMax*(1 - u / densityJam);
}

double ve(double u)
{
	return vMax*(1 / (1 + exp((u / densityJam - 0.25) / 0.06)) - 3.72e-6);
}

double intOmega(double u)
{
	return (1 - ve(u) / vMax)*densityJam;
}

//流量密度函数
double flow(double u)
{
	return u*V(u);
}

//Riemann需求函数
double demand(double u)
{
	return u > densityCr ? flowCr : flow(u);
}

//Riemann供给函数
double supply(double u)
{
	return u < densityCr ? flowCr : flow(u);
}

//Godunov数值流通量
double FGodunov(double ul, double ur)
{
	if (ul < densityCr)
	{
		if (ur <= densityCr)
			return flow(ul);
		else
		{
			if ((flow(ur) - flow(ul))*(ur - ul) >= 0)
				return flow(ul);
			else
				return flow(ur);
		}
	}
	else
	{
		if (ur>densityCr)
			return flow(ur);
		else
			return flowCr;
	}
}

double veLaneN(double rho, int n)
{
	return vMax*(1 / (1 + exp((rho / (n*densityJam) - 0.25) / 0.06)) - 3.72e-6);
}

double VLaneN(double rho, int n)
{
	return vMax*(1 - rho / (n*densityJam));
}

double flowLaneN(double rho, int n)
{
	return rho*VLaneN(rho, n);
}

double FGLaneN(double rho1, double rho2, int n)
{
	if (rho1<densityJam*n / 2)
	{
		if (rho2<densityJam*n / 2)
			return flowLaneN(rho1, n);
		else
		{
			if ((flowLaneN(rho2, n) - flowLaneN(rho1, n))*(rho2 - rho1) >= 0)
				return flowLaneN(rho1, n);
			else
				return flowLaneN(rho2, n);
		}
	}
	else
	{
		if (rho2>densityJam*n / 2)
			return flowLaneN(rho2, n);
		else
			return flowLaneN(n*densityJam / 2, n);
	}
}

double demandN(double rho, int n)
{
	if (rho<n*densityJam / 2)
		return flowLaneN(rho, n);
	else
		return flowLaneN(n*densityJam / 2, n);
}

double supplyN(double rho, int n)
{
	if (rho<n*densityJam / 2)
		return flowLaneN(n*densityJam / 2, n);
	else
		return flowLaneN(rho, n);
}


//向量逆时针旋转alpha弧度
void trans(double &x, double &y, double alpha)
{
	x = x*cos(alpha) - y*sin(alpha);
	y = x*sin(alpha) + y*cos(alpha);
}

//全局初始化函数
void ini(Element* element)
{
	//特殊车道数赋值
	for (int i = 866; i < 881; i++)
	{
		element[i].laneN = 8;
	}

	//密度初值
	for (int i = 0; i < N; i++)
	{
		if (i < N1)
			element[i].initialization(0.2*densityJam*element[i].laneN);//Road1初值
		else if (i < N1 + N2)
			element[i].initialization(0.2*densityJam*element[i].laneN);//Road2初值
		else
			element[i].initialization(0.2*densityJam*element[i].laneN);//Road3初值

		element[i].setValueXY(Det_x / 2 + Det_x*i, 0);

		//		if (i >= 200 && i < 400)//微小正弦扰动
		//			element[i].density += 0.1*densityJam*(sin(2 * PI / 4000 * (element[i].getValueX() - element[200].getValueX())));
	}

	//求网格位置
	double x2 = Length12*cos(alpha1) + Length11;//拐角2坐标x
	double y2 = Length12*sin(alpha1);//拐角2坐标y

	double x3 = Length13*cos(alpha2) + x2;//拐角3坐标x
	double y3 = Length13*sin(alpha2) + y2;//拐角3坐标y

	double x4 = Length21*cos(alpha3) + x3;//拐角4坐标x
	double y4 = Length21*sin(alpha3) + y3;//拐角4坐标y

	double x5 = Length22*cos(alpha4) + x4;//拐角5坐标x
	double y5 = Length22*sin(alpha4) + y4;//拐角5坐标y

	double x6 = Length31*cos(alpha5) + x3;//拐角6坐标x
	double y6 = Length31*sin(alpha5) + y3;//拐角6坐标y

	double x7 = Length32*cos(alpha6) + x6;//拐角7坐标x
	double y7 = Length32*sin(alpha6) + y6;//拐角7坐标y

	double x0;
	for (int i = N11; i < N; i++)//旋转改变网格中心位置
	{
		if (i < N11 + N12)
		{
			x0 = element[i].getValueX() - Length11;
			element[i].setValueXY(x0*cos(alpha1) + Length11, x0*sin(alpha1) + 0);
		}
		else if (i < N1)
		{
			x0 = element[i].getValueX() - (Length11 + Length12);
			element[i].setValueXY(x0*cos(alpha2) + x2, x0*sin(alpha2) + y2);
		}
		else if (i < N1 + N21)//第一个分叉
		{
			x0 = element[i].getValueX() - Length1;
			element[i].setValueXY(x0*cos(alpha3) + x3, x0*sin(alpha3) + y3);
		}
		else if (i < N1 + N2)
		{
			x0 = element[i].getValueX() - (Length1 + Length21);
			element[i].setValueXY(x0*cos(alpha4) + x4, x0*sin(alpha4) + y4);
		}
		else if (i < N1 + N2 + N31)
		{
			x0 = element[i].getValueX() - (Length1 + Length2);
			element[i].setValueXY(x0*cos(alpha5) + x3, x0*sin(alpha5) + y3);
		}
		else if (i < N1 + N2 + N31 + N32)
		{
			x0 = element[i].getValueX() - (Length1 + Length2 + Length31);
			element[i].setValueXY(x0*cos(alpha6) + x6, x0*sin(alpha6) + y6);
		}
		else if (i < N)
		{
			x0 = element[i].getValueX() - (Length1 + Length2 + Length31 + Length32);
			element[i].setValueXY(x0*cos(alpha7) + x7, x0*sin(alpha7) + y7);
		}
	}
	//特殊flag赋值
	element[0].flag = 0;//进口
	element[N1 + N2 - 1].flag = 1;//出口
	element[N - 1].flag = 1;//出口
	element[N1 - 1].flag = 2;//分岔口前网格
	element[N1].flag = 3;//分岔口后网格
	element[N1 + N2].flag = 3;//分岔口后网格

}

//更新数值流通量函数
void updateFlux(Element *element)
{
	for (int i = 0; i < N; i++)
	{
		if (element[i].flag == 0)//进口边界网格、要改成固定边界条件
		{
			element[i].flux = FGLaneN(element[i].density, element[element[i].EOut[0]].density, element[i].laneN);//Godunov数值流通量
			element[element[i].EOut[0]].flowIn = element[i].flux;
			element[i].flowIn = FGLaneN(element[i].density, element[i].density, element[i].laneN);//自然边界条件的Godunov数值流通量
		}
		else if (element[i].flag == 1)//出口边界网格
		{
			element[i].flux = FGLaneN(element[i].density, element[i].density, element[i].laneN);//自然边界条件的Godunov数值流通量
		}
		else if (element[i].flag == 2)//分岔口网格
		{
			double d = demandN(element[i].density, element[i].laneN);
			double s0 = supplyN(element[element[i].EOut[0]].density, element[element[i].EOut[0]].laneN);
			double s1 = supplyN(element[element[i].EOut[1]].density, element[element[i].EOut[1]].laneN);
			element[i].flux = std::min(std::min(d, s0 / RateFlowDiverge), s1 / (1 - RateFlowDiverge));
			element[element[i].EOut[0]].flowIn = RateFlowDiverge*element[i].flux;
			element[element[i].EOut[1]].flowIn = (1 - RateFlowDiverge)*element[i].flux;

		}
		else//一般网格
		{
			double d = demandN(element[i].density, element[i].laneN);
			double s = supplyN(element[element[i].EOut[0]].density, element[element[i].EOut[0]].laneN);
			element[i].flux = std::min(d, s);//FGLaneN(element[i].density, element[element[i].EOut[0]].density, element[i].laneN);//Godunov数值流通量
			element[element[i].EOut[0]].flowIn = element[i].flux;
		}
	}
}

//上下游网格更新函数
void updateEIO(Element* element)
{
	for (int i = 0; i < N; i++)
	{
		if (i == N1 - 1)
		{
			element[i].EOut[0] = N1;
			element[i].EOut[1] = N1 + N2;
		}
		else if (i == N1 + N2 - 1 || i == N - 1)
		{
			element[i].EOut[0] = -1;
			element[i].EOut[1] = -1;
		}
		else
		{
			element[i].EOut[0] = i + 1;
			element[i].EOut[1] = -1;
		}
	}
}

//自定义时间函数
double cal_g(double t)
{
	//
	double G1 = 0.8, G2 = 0.6;
	double t1 = 70, t2 = 90, t3 = 100, t4 = 120, t5 = 150;
	if (t <= 0)
		return 0;
	else if (t <= t1)
		return G1*t / t1;
	else if (t <= t2)
		return G1;
	else if (t <= t3)
		return (G2 - G1)*(t - t2) / (t3 - t2) + G1;
	else if (t <= t4)
		return G2;
	else if (t <= t5)
		return G2*(t - t5) / (t4 - t5);
	else
		return 0;
}

//流入流出更改函数，模拟匝道io
void updateST(Element* element, double t)//改成右端项,移项合并而来
{
	for (int i = 0; i < NI; i++)
	{
		int iNE = int(NIx[i] / Det_x);
		element[iNE].rampIn = cal_g(t)*0.2*densityJam*element[i].laneN*(1 - element[i].density / (densityJam*element[i].laneN));
	}
	for (int i = 0; i < NO; i++)
	{
		int oNE = int(NOx[i] / Det_x);
		element[oNE].rampOut = 0.1*element[oNE].flux / Det_x;//+= cal_g(t-15)*0.2*densityJam*element[i].density / densityJam*Det_x*Det_x;
	}
}

//网格密度更新函数
double det_t;
void updateElement(Element *element, double t)
{
	updateFlux(element);//更新流量
	updateST(element, t);//更新源项

	det_t = 1;//CFL_Det_t(element);//根据CFL条件求最大允许时间步
	if (!det_t)
	{
		printf("-----------delta_t = 0，程序终止-----------\n");
		system("pause");
	}
	if (t - 5 * floor(t / 5.0) <= 1)//时间函数，使收费站处流量强行为0
	{
		printf("%f\n", t);
		element[874].flux = 0;
		element[875].flowIn = 0;
	}
	for (int i = 1; i < N; i++)
	{
		element[i].density = element[i].density + det_t / Det_x*(element[i].flowIn - element[i].flux) + det_t*(element[i].rampIn - element[i].rampOut);
	}
}

//非齐次边界条件
double CFL_Det_t_ST(Element* element)
{
	double dt = 100;// CFL_Det_t(element);
	for (int i = 0; i < N; i++)
	{
		double density1 = 0, density2 = element[i].density, density3 = 0;
		if (element[i].flag == 0)//进口
		{
			density1 = element[i].density;
			density3 = element[element[i].EOut[0]].density;
		}
		else if (element[i].flag == 1)//出口
		{
			density1 = element[i - 1].density;
			density3 = element[i].density;
		}
		else if (element[i].flag == 2)//交叉口前网格
		{
			density1 = element[i - 1].density;
			density3 = (densityJam + sqrt(densityJam*densityJam - 4 * element[i].flux*densityJam / vMax)) / 2;
		}
		else if (element[i].flag == 3)//交叉口后网格
		{
			density1 = (densityJam - sqrt(densityJam*densityJam - 4 * element[i].flowIn*densityJam / vMax)) / 2;
			density3 = element[element[i].EOut[0]].density;
		}
		else//其他正常网格
		{
			density1 = element[i - 1].density;
			density3 = element[element[i].EOut[0]].density;
		}
		double Min = std::min(std::min(density1, density2), density3);
		double Max = std::max(std::max(density1, density2), density3);

		if (element[i].rampIn<element[i].rampOut)
		{
			if (dt > fabs(-Min / (element[i].rampIn - element[i].rampOut)))
				dt = fabs(-Min / (element[i].rampIn - element[i].rampOut));
		}
		else if (element[i].rampIn > element[i].rampOut)
		{
			if (dt > fabs((densityJam - Max) / (element[i].rampIn - element[i].rampOut)))
				dt = fabs((densityJam - Max) / (element[i].rampIn - element[i].rampOut));
		}
		else
		{
			dt = 0.1;
		}
	}
	return dt;
}
//齐次方程CFL条件
double CFL_Det_t(Element* element)
{
	double maxDf = abs(vMax*(1 - 2 * element[0].density / (densityJam*element[0].laneN)));
	for (int i = 1; i < N; i++)
		if (abs(vMax*(1 - 2 * element[i].density / (densityJam*element[i].laneN)))>maxDf)
			maxDf = abs(vMax*(1 - 2 * element[i].density / (densityJam*element[i].laneN)));
	if (maxDf)
		return Det_x / maxDf;
	else
		return 0.1;
}

//yzd自拟定CFL条件
double CFL_Det_t_yzd(Element* element)
{
	double maxDf = abs((element[0].flux - element[0].flowIn) / element[0].density);
	for (int i = 1; i < N; i++)
		if (abs((element[i].flux - element[i].flowIn) / element[i].density)>maxDf)
			maxDf = abs((element[i].flux - element[i].flowIn) / element[i].density);
	if (maxDf)
		return Det_x / maxDf - 0.0001;
	else
		return 1;//这里到底应该返回什么？
}
