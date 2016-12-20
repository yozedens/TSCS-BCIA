#include"pltData.h"

#ifndef   SETVALUE_H

#include"SetValue.h"
#endif

#include<fstream>
#include<iostream>

using namespace std;

void coutChangedData(Element* element)
{

	for (int i = 0; i < N; i++)
	{
		if (abs(element[i].density - 0.05)>0.0001)
			std::cout << i << ":" << element[i].density << "  ";
	}
	std::cout << std::endl;
}

void pltTSDensity2FileforUnity(Element* element, int time)//���һ���ļ�����Unityʹ��
{
	char tempChar[15];
	sprintf_s(tempChar, "%d", time);
	string timeString = tempChar;

	ofstream fileDensity(".\\out\\Density" + timeString + ".txt", ios_base::app);// ios_base::app��ʾ���������ļ������ļ�β׷�� ���û���򴴽�

	for (int i = 0; i < N - 1; i++)
	{
		fileDensity << element[i].density / element[i].laneN << ",";
	}
	fileDensity << element[N - 1].density / element[N - 1].laneN << endl;
	fileDensity.close();

}

void pltTSDensity2FileforTecplot(Element* element, int time, double t)//���һ���ļ�����������ƽ��
{
	char tempChar[15];
	sprintf_s(tempChar, "%d", time);
	string timeString = tempChar;

	ofstream fileDensity1(".\\out\\densityRoad"+ timeString+".plt", ios_base::app);// ios_base::app��ʾ���������ļ������ļ�β׷�� ���û���򴴽�

	fileDensity1 << "          zone t=    \"" << t << "\"" << endl;
	for (int i = 0; i < N; i++)
	{
		fileDensity1 << i << " " << element[i].density / element[i].laneN << endl;
	}
	fileDensity1.close();
}

void pltGridPosition2File(Element* element)
{
	ofstream fileGrid(".\\out\\Grid.plt");

	for (int i = 0; i < N; i++)
	{
		fileGrid << element[i].getValueX() << " " << element[i].getValueY() << endl;
	}
	//fileGrid << element[N - 1].getValueX() << "," << element[N - 1].getValueY() << "," << 0.0 << endl;

	fileGrid.close();
}

void pltGridPosition2FileforUnity(Element* element)
{
	ofstream fileGrid(".\\out\\GridforUnity.plt");

	for (int i = 0; i < N - 1; i++)
	{
		fileGrid << element[i].getValueX() << "," << element[i].getValueY() << ",";
	}
	fileGrid << element[N - 1].getValueX() << "," << element[N - 1].getValueY() << endl;

	fileGrid.close();
}