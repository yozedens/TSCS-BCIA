#include"SetValue.h"
#include<algorithm>


double veLaneN(double rho, int n)
{
	return vMax*(1 / (1 + exp((rho / (n*densityJam) - 0.25) / 0.06)) - 3.72e-6);
}

double VLaneN(double rho, int n)
{
	return vMax*(1 - rho / (n*densityJam));
}

double iniOmega(double u, int n)
{
	return (1 - veLaneN(u, n) / vMax)*densityJam*n;
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


//������ʱ����תalpha����
void trans(double &x, double &y, double alpha)
{
	x = x*cos(alpha) - y*sin(alpha);
	y = x*sin(alpha) + y*cos(alpha);
}

//ȫ�ֳ�ʼ������
void ini(Element* element, int UorD)
{
	//���⳵������ֵ
	for (int i = 866; i < 881; i++)
	{
		if (UorD)
			element[i].laneN = 5;
		else
			element[i].laneN = 8;
	}

	//�ܶȳ�ֵ
	for (int i = 0; i < N; i++)
	{
		if (i < N1)
			element[i].initialization(0.2*densityJam*element[i].laneN);//Road1��ֵ
		else if (i < N)
			element[i].initialization(0.2*densityJam*element[i].laneN);//Road2��ֵ

		element[i].omega = iniOmega(element[i].density, element[i].laneN);

		element[i].setValueXY(Det_x / 2 + Det_x*i, 0);

		//		if (i >= 200 && i < 400)//΢С�����Ŷ�
		//			element[i].density += 0.1*densityJam*(sin(2 * PI / 4000 * (element[i].getValueX() - element[200].getValueX())));
	}

	//������λ��
	double x2 = Length12*cos(alpha1) + Length11;//�ս�2����x
	double y2 = Length12*sin(alpha1);//�ս�2����y

	double x3 = Length13*cos(alpha2) + x2;//�ս�3����x
	double y3 = Length13*sin(alpha2) + y2;//�ս�3����y

	double x4 = Length21*cos(alpha3) + x3;//�ս�4����x
	double y4 = Length21*sin(alpha3) + y3;//�ս�4����y

	double x5 = Length22*cos(alpha4) + x4;//�ս�5����x
	double y5 = Length22*sin(alpha4) + y4;//�ս�5����y


	double x0;
	for (int i = N11; i < N; i++)//��ת�ı���������λ��
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
		else if (i < N1 + N21)
		{
			x0 = element[i].getValueX() - Length1;
			element[i].setValueXY(x0*cos(alpha3) + x3, x0*sin(alpha3) + y3);
		}
		else if (i < N)
		{
			x0 = element[i].getValueX() - (Length1 + Length21);
			element[i].setValueXY(x0*cos(alpha4) + x4, x0*sin(alpha4) + y4);
		}
	}
	//����flag��ֵ��ȱʡֵΪ10
	if (UorD)
	{
		element[N - 1].flag = 0;//����
		element[0].flag = 1;//����
	}
	else
	{
		element[0].flag = 0;//����
		element[N - 1].flag = 1;//����
	}

}

//������ֵ��ͨ������
void updateFlux(Element *element)
{
	for (int i = 0; i < N; i++)
	{
		if (element[i].flag == 0)//���ڱ߽�����Ҫ�ĳɹ̶��߽�����
		{
			element[i].flux = FGLaneN(element[i].omega, element[element[i].EOut[0]].omega, element[i].laneN);//Godunov��ֵ��ͨ��
			element[element[i].EOut[0]].flowIn = element[i].flux;
			element[i].flowIn = FGLaneN(iniOmega(boundryDensity*element[i].laneN, element[i].laneN), element[i].omega, element[i].laneN);//�̶��߽�������Godunov��ֵ��ͨ��
/*			double d = demandN(iniOmega(boundryDensity*element[i].laneN, element[i].laneN), element[i].laneN);
			double s = supplyN(element[i].omega, element[i].laneN);
			element[i].flowIn = std::min(d, s);
			
			d = demandN(element[i].omega, element[i].laneN);
			s = supplyN(element[element[i].EOut[0]].omega, element[element[i].EOut[0]].laneN);
			element[i].flux = std::min(d, s);
			element[element[i].EOut[0]].flowIn = element[i].flux;
*/
		}
		else if (element[i].flag == 1)//���ڱ߽�����
		{
			element[i].flux = FGLaneN(element[i].omega, element[i].omega, element[i].laneN);//��Ȼ�߽�������Godunov��ֵ��ͨ��
		}
		else if (element[i].flag == 2)//�ֲ������
		{
			double d = demandN(element[i].omega, element[i].laneN);
			double s0 = supplyN(element[element[i].EOut[0]].omega, element[element[i].EOut[0]].laneN);
			double s1 = supplyN(element[element[i].EOut[1]].omega, element[element[i].EOut[1]].laneN);
			element[i].flux = std::min(std::min(d, s0 / RateFlowDiverge), s1 / (1 - RateFlowDiverge));
			element[element[i].EOut[0]].flowIn = RateFlowDiverge*element[i].flux;
			element[element[i].EOut[1]].flowIn = (1 - RateFlowDiverge)*element[i].flux;

		}
		else//һ������
		{
			double d = demandN(element[i].omega, element[i].laneN);
			double s = supplyN(element[element[i].EOut[0]].omega, element[element[i].EOut[0]].laneN);
			element[i].flux = std::min(d, s);
			element[element[i].EOut[0]].flowIn = element[i].flux;
		}
	}
}

//����������º���
void updateEIO(Element* element, int UorD)
{
	for (int i = 0; i < N; i++)
	{
		if (element[i].flag==1)
		{
			element[i].EOut[0] = -1;
			element[i].EOut[1] = -1;
		}
		else
		{
			if (UorD)//Ϊ�Ż������ٶȣ��˴��жϿɸ�Ϊһ�����жϣ��������
			{
				element[i].EOut[0] = i - 1;	
			}
			else
			{
				element[i].EOut[0] = i + 1;
			}
			element[i].EOut[1] = -1;
		}
	}
}

//�Զ���ʱ�亯��
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

//�����������ĺ�����ģ���ѵ�io
void updateST(Element* element, double t, int UorD)//�ĳ��Ҷ���,����ϲ�����
{
	for (int i = 0; i < NI; i++)
	{
		int iNE = int(NIx[i] / Det_x)*(1 - UorD) + int(NIDx[i] / Det_x)*UorD;
		element[iNE].rampIn = cal_g(t)*0.1*densityJam*element[i].laneN*(1 - element[i].density / (densityJam*element[i].laneN));
	}
	for (int i = 0; i < NO; i++)
	{
		int oNE = int(NOx[i] / Det_x)*(1 - UorD) + int(NODx[i] / Det_x)*UorD;
		element[oNE].rampOut = 0.1*element[oNE].flux / Det_x;//+= cal_g(t-15)*0.2*densityJam*element[i].density / densityJam*Det_x*Det_x;
	}
}

//�����ܶȸ��º���
double det_t;
void updateElement(Element *element, double t, int UorD)
{
	updateFlux(element);//��������
	updateST(element, t, UorD);//����Դ��

	det_t = 0.8;//CFL_Det_t(element);//����CFL�������������ʱ�䲽
	if (!det_t)
	{
		printf("-----------delta_t = 0��������ֹ-----------\n");
		system("pause");
	}
	for (int i = 0; i < N; i++)//����omega
	{
		element[i].omega = element[i].omega + det_t / Det_x*(element[i].flowIn - element[i].flux)
			+ det_t*((VLaneN(element[i].omega, element[i].laneN) - veLaneN(element[i].density, element[i].laneN)) / beta);
	}

	if (t - 5 * floor(t / 5.0) <= 1)//ʱ�亯����ʹ�շ�վ������ǿ��Ϊ0
	{
		Fb(element, UorD);
		//		printf("%f\n", t);		
	}

	if (UorD)//����density
	{
		for (int i = 0; i < N - 1; i++)
		{
			element[i].density = element[i].density + det_t / Det_x*(element[i + 1].density / element[i + 1].omega*element[i].flowIn
				- element[i].density / element[i].omega*element[i].flux) + det_t*(element[i].rampIn - element[i].rampOut);
		}
		element[N - 1].density = element[N - 1].density + det_t / Det_x*(boundryDensity*element[N - 1].laneN / iniOmega(boundryDensity*element[N - 1].laneN, element[N - 1].laneN)*element[N - 1].flowIn
			- element[N - 1].density / element[N - 1].omega*element[N - 1].flux) + det_t*(element[N - 1].rampIn - element[N - 1].rampOut);
	}
	else
	{
		element[0].density = element[0].density + det_t / Det_x*(boundryDensity*element[0].laneN / iniOmega(boundryDensity*element[0].laneN, element[0].laneN)*element[0].flowIn
			- element[0].density / element[0].omega*element[0].flux) + det_t*(element[0].rampIn - element[0].rampOut);
		for (int i = 1; i < N; i++)
		{
			element[i].density = element[i].density + det_t / Det_x*(element[i - 1].density / element[i - 1].omega*element[i].flowIn
				- element[i].density / element[i].omega*element[i].flux) + det_t*(element[i].rampIn - element[i].rampOut);
		}
	}
}

void Fb(Element* element, int UorD)//������������
{
	if (UorD)
	{
		element[875].flux = 0;
		element[874].flowIn = 0;
	}
	else
	{
		element[874].flux = 0;
		element[875].flowIn = 0;
	}
}
//����α߽�����
double CFL_Det_t_ST(Element* element)
{
	double dt = 100;// CFL_Det_t(element);
	for (int i = 0; i < N; i++)
	{
		double density1 = 0, density2 = element[i].density, density3 = 0;
		if (element[i].flag == 0)//����
		{
			density1 = element[i].density;
			density3 = element[element[i].EOut[0]].density;
		}
		else if (element[i].flag == 1)//����
		{
			density1 = element[i - 1].density;
			density3 = element[i].density;
		}
		else if (element[i].flag == 2)//�����ǰ����
		{
			density1 = element[i - 1].density;
			density3 = (densityJam + sqrt(densityJam*densityJam - 4 * element[i].flux*densityJam / vMax)) / 2;
		}
		else if (element[i].flag == 3)//����ں�����
		{
			density1 = (densityJam - sqrt(densityJam*densityJam - 4 * element[i].flowIn*densityJam / vMax)) / 2;
			density3 = element[element[i].EOut[0]].density;
		}
		else//������������
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
//��η���CFL����
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

//yzd���ⶨCFL����
double CFL_Det_t_yzd(Element* element)
{
	double maxDf = abs((element[0].flux - element[0].flowIn) / element[0].density);
	for (int i = 1; i < N; i++)
		if (abs((element[i].flux - element[i].flowIn) / element[i].density)>maxDf)
			maxDf = abs((element[i].flux - element[i].flowIn) / element[i].density);
	if (maxDf)
		return Det_x / maxDf - 0.0001;
	else
		return 1;//���ﵽ��Ӧ�÷���ʲô��
}
