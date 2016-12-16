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

void pltTSDensity2File(Element* element, int time)//输出一个文件，供Unity使用
{
	char tempChar[15];
	sprintf_s(tempChar, "%d", time);
	string timeString = tempChar;

	ofstream fileDensity(".\\out\\Density" + timeString + ".txt", ios_base::app);// ios_base::app表示如果有这个文件则在文件尾追加 如果没有则创建

	for (int i = 0; i < N - 1; i++)
	{
		fileDensity << element[i].density / element[i].laneN << ",";
	}
	fileDensity << element[N - 1].density / element[N - 1].laneN << endl;
	fileDensity.close();

}

void pltTSDensity2File2(Element* element, double t)
{
	char tempChar[15];
	sprintf_s(tempChar, "%f", t);//将double型转化为string型的关键步骤，%f表示double型
	string tempString = tempChar;
	string outfileName1 = ".\\out\\pltRho";
	outfileName1 = outfileName1 + tempString + "_1.plt";
	ofstream fileDensity1(outfileName1);

	for (int i = 0; i < N1; i++)
	{
		fileDensity1 << i << " " << element[i].density / element[i].laneN << endl;
	}

	fileDensity1.close();

	string outfileName2 = ".\\out\\pltRho";
	outfileName2 = outfileName2 + tempString + "_2.plt";
	ofstream fileDensity2(outfileName2);

	for (int i = N1; i < N1 + N2; i++)
	{
		fileDensity2 << i - N1 << " " << element[i].density / element[i].laneN << endl;
	}

	fileDensity2.close();

	string outfileName3 = ".\\out\\pltRho";
	outfileName3 = outfileName3 + tempString + "_3.plt";
	ofstream fileDensity3(outfileName3);

	for (int i = N1 + N2; i < N; i++)
	{
		fileDensity3 << i - N1 - N2 << " " << element[i].density / element[i].laneN << endl;
	}

	fileDensity3.close();

}

void pltTSDensity2File3(Element* element, double t)//每条路输出一个文件，已做车道平均
{

	ofstream fileDensity1(".\\out\\densityRoad1.plt", ios_base::app);// ios_base::app表示如果有这个文件则在文件尾追加 如果没有则创建

	fileDensity1 << "          zone t=    \"" << t << "\"" << endl;
	for (int i = 0; i < N1; i++)
	{
		fileDensity1 << i << " " << element[i].density / element[i].laneN << endl;
	}
	fileDensity1.close();


	ofstream fileDensity2(".\\out\\densityRoad2.plt", ios_base::app);// ios_base::app表示如果有这个文件则在文件尾追加 如果没有则创建
	fileDensity2 << "          zone t=    \"" << t << "\"" << endl;
	for (int i = N1; i < N1 + N2; i++)
	{
		fileDensity2 << i - N1 << " " << element[i].density / element[i].laneN << endl;
	}
	fileDensity2.close();

	ofstream fileDensity3(".\\out\\densityRoad3.plt", ios_base::app);// ios_base::app表示如果有这个文件则在文件尾追加 如果没有则创建
	fileDensity3 << "          zone t=    \"" << t << "\"" << endl;
	for (int i = N1 + N2; i < N; i++)
	{
		fileDensity3 << i - N1 - N2 << " " << element[i].density / element[i].laneN << endl;
	}
	fileDensity3.close();
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