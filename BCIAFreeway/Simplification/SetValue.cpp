#include"SetValue.h"
#include<algorithm>

double ve(double u)//���ܹ�ϵ����
{
	return vMax*(1 - u / densityJam);
}

double flow(double u)//��������
{
	return u*ve(u);
}

double demand(double u)
{
	return u > densityCr ? flowCr : flow(u);
}

double supply(double u)
{
	return u < densityCr ? flowCr : flow(u);
}

double FGodunov(double ul, double ur)//Godunov��ֵ��ͨ��
{
	if (ul < densityCr)
	{
		if (ur <= densityCr)
			return flow(ul);
		else
		{
			if ((flow(ur)-flow(ul))*(ur-ul) >= 0)
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

void trans(double &x, double &y, double alpha)//������ת
{
	x = x*cos(alpha) - y*sin(alpha);
	y = x*sin(alpha) + y*cos(alpha);
}

void ini(Element* element)
{
	for (int i = 0; i < N; i++)
	{
		if (i < N1)
			element[i].initialization(0.6*densityJam);
		else
			element[i].initialization(0.4*densityJam);

		element[i].setValueXY(Det_x / 2 + Det_x*i, 0);

		if (i >= 200 && i < 400)//΢С�����Ŷ�
			element[i].density += 0.2*densityJam*(sin(2 * PI / 4000 * (element[i].getValueX() - element[200].getValueX())));
	}
	

	double x2 = Length12*cos(alpha1) + Length11;//�ս�2����x
	double y2 = Length12*sin(alpha1) ;//�ս�2����y

	double x3 = Length13*cos(alpha2) + x2;//�ս�3����x
	double y3 = Length13*sin(alpha2) + y2;//�ս�3����y

	double x4 = Length21*cos(alpha3) + x3;//�ս�4����x
	double y4 = Length21*sin(alpha3) + y3;//�ս�4����y

	double x5 = Length22*cos(alpha4) + x4;//�ս�5����x
	double y5 = Length22*sin(alpha4) + y4;//�ս�5����y

	double x6 = Length31*cos(alpha5) + x3;//�ս�6����x
	double y6 = Length31*sin(alpha5) + y3;//�ս�6����y

	double x7 = Length32*cos(alpha6) + x6;//�ս�7����x
	double y7 = Length32*sin(alpha6) + y6;//�ս�7����y

/*	for (int i = N11; i < N; i++)//��ת�ı���������λ��
	{
		x0 = element[i].getValueX() - x1;
		y0 = element[i].getValueY() - y1;
		trans(x0, y0, alpha1);
		element[i].setValueXY(x0 + x1, y0 + y1);
		if (i >= N11 + N12)
		{
			x0 = element[i].getValueX() - x2;
			y0 = element[i].getValueY() - y2;
			trans(x0, y0, alpha2);
			element[i].setValueXY(x0 + x2, y0 + y2);
			if (i >= N1&&i<N1+N2)//��һ���ֲ�
			{
				x0 = element[i].getValueX() - x3;
				y0 = element[i].getValueY() - y3;
				trans(x0, y0, alpha3);
				element[i].setValueXY(x0 + x3, y0 + y3);
				if (i >= N1 + N21)
				{
					x0 = element[i].getValueX() - x4;
					y0 = element[i].getValueY() - y4;
					trans(x0, y0, alpha4);
					element[i].setValueXY(x0 + x4, y0 + y4);
				}
			}
			else//�ڶ����ֲ�
			{
				x0 = element[i].getValueX() - x5;
				y0 = element[i].getValueY() - y5;
				trans(x0, y0, alpha5);
				element[i].setValueXY(x0 + x5, y0 + y5);
				if (i >= N1 + N2 + N31)
				{
					x0 = element[i].getValueX() - x6;
					y0 = element[i].getValueY() - y6;
					trans(x0, y0, alpha6);
					element[i].setValueXY(x0 + x6, y0 + y6);
					if (i >= N1 + N2 + N32)
					{
						x0 = element[i].getValueX() - x7;
						y0 = element[i].getValueY() - y7;
						trans(x0, y0, alpha7);
						element[i].setValueXY(x0 + x7, y0 + y7);
					}
				}
			}
		}
	}
	*/
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
		else if (i < N1 + N21)//��һ���ֲ�
		{
			x0 = element[i].getValueX() - Length1;
			element[i].setValueXY(x0*cos(alpha3) + x3, x0*sin(alpha3) + y3);
		}
		else if (i < N1 + N2)
		{
			x0 = element[i].getValueX() - (Length1 + Length21);
			element[i].setValueXY(x0*cos(alpha4) + x4, x0*sin( alpha4) + y4);
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

	element[0].flag = 0;//����
	element[N1 + N2 - 1].flag = 1;//����
	element[N - 1].flag = 1;//����
	element[N1 - 1].flag = 2;//�ֲ��
}

void updateFlux(Element *element)//������ֵ��ͨ��
{
	for (int i = 0; i < N; i++)
	{
		if (element[i].flag == 0)//���ڱ߽�����
		{
			element[i].flux = FGodunov(element[i].density, element[element[i].EOut[0]].density);//Godunov��ֵ��ͨ��
			element[element[i].EOut[0]].flowIn = element[i].flux;
			element[i].flowIn = FGodunov(element[i].density, element[i].density);//��Ȼ�߽�������Godunov��ֵ��ͨ��
		}
		else if (element[i].flag == 1)//���ڱ߽�����
		{
			element[i].flux = FGodunov(element[i].density, element[i].density);//��Ȼ�߽�������Godunov��ֵ��ͨ��
		}
		else if (element[i].flag == 2)//�ֲ������
		{
			double d = demand(element[i].density);
			double s0 = supply(element[element[i].EOut[0]].density);
			double s1 = supply(element[element[i].EOut[1]].density);
			element[i].flux = std::min(std::min(d,s0/RateFlowDiverge), s1/(1-RateFlowDiverge));
//			double min = d;
//			if (s0 / RateFlowDiverge<min)
//				min = s0 / RateFlowDiverge;
//			if (s1 / (1.0 - RateFlowDiverge)<min)
//				min = s1 / (1.0 - RateFlowDiverge);
//			element[i].flux = min;

			element[element[i].EOut[0]].flowIn = RateFlowDiverge*element[i].flux;
			element[element[i].EOut[1]].flowIn = (1 - RateFlowDiverge)*element[i].flux;

		}
		else//һ������
		{
			element[i].flux = FGodunov(element[i].density, element[element[i].EOut[0]].density);//Godunov��ֵ��ͨ��
			element[element[i].EOut[0]].flowIn = element[i].flux;
		}
	}
}

void updateEIO(Element* element)
{
	for (int i = 0; i < N; i++)
	{
		if (i == N1-1)
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

void updateST(Element* element)
{
	for (int i = 0; i < NI; i++)
	{
		int iNE = int(NIx[i]/Det_x);
		element[iNE].flowIn *= 1.1;
	}
	for (int i = 0; i < NO; i++)
	{
		int oNE = int(NOx[i] / Det_x);
		element[oNE].flux *= 1.1;
	}
}

double det_t;

void updateElement(Element *element)
{
	updateFlux(element);//��������
	updateST(element);//����Դ��
	det_t = 0.01;// CFL_Det_t(element);//����CFL�������������ʱ�䲽
	if (!det_t)
	{
		printf("-----------delta_t = 0��������ֹ-----------\n");
		system("pause");
	}
	for (int i = 0; i < N; i++)
	{
		if (true)
		{
			element[i].density = element[i].density + det_t / Det_x*(element[i].flowIn - element[i].flux);
		}
	}
}

double CFL_Det_t(Element* element)
{
	double maxDf = abs((element[0].flux - element[0].flowIn) / element[0].density);
	for (int i = 1; i < N; i++)
		if (abs((element[i].flux - element[i].flowIn) / element[i].density)>maxDf)
			maxDf = abs((element[i].flux - element[i].flowIn) / element[i].density);
	if (maxDf)
		return Det_x / maxDf-0.0001;
	else
		return 1;//���ﵽ��Ӧ�÷���ʲô��
}
double CFL_Det_t1(Element* element)
{
	double maxDf = abs(vMax*(1 - 2 * element[0].density / densityJam));
	for (int i = 1; i < N; i++)
		if (abs(vMax*(1 - 2 * element[i].density / densityJam))>maxDf)
			maxDf = abs(vMax*(1 - 2 * element[i].density / densityJam));
	if (maxDf)
		return Det_x / maxDf;
	else
		return 0;
}
