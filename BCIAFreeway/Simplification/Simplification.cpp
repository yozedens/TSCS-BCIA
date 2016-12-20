#include<iostream>
#include<fstream>
#include"NormalTFClass.h"
#include"SetValue.h"
#include"pltData.h"
#include<time.h>
using namespace std;

int main()
{
	Element element[N];
	ini(element);

	double t = 0;

	/*	ofstream fileGt(".\\out\\gt.plt");
	for (double tt = 0; tt < T; tt += 1)
	fileGt << tt << " " << cal_g(tt) << endl;
	fileGt.close();
	system("pause");
	*/
	int timeInt = time(0);
	pltTSDensity2FileforUnity(element, timeInt);
	pltTSDensity2FileforTecplot(element, timeInt, 0);
//	pltGridPosition2File(element);
//	pltGridPosition2FileforUnity(element);
	updateEIO(element);

	while (t < T)
	{
		t += det_t;
		//		if(count%10==0)
		std::cout << "t= " << t << ":" << element[N1 - 1].flux << "  " << element[N1].flowIn << " " << element[N1 + N2].flowIn << std::endl;
		updateElement(element, t);
		//		coutChangedData(element);
		//		if (count % 10 == 0)
		{
			pltTSDensity2FileforUnity(element, timeInt);
			pltTSDensity2FileforTecplot(element, timeInt, t);
		}
	}

	Element elementDown[N];
	timeInt = time(0);
	t = 0;
	ini(elementDown, 1);
	pltTSDensity2FileforUnity(elementDown, timeInt);
	pltTSDensity2FileforTecplot(elementDown, timeInt, 0);
	updateEIO(elementDown, 1);

	while (t < T)
	{
		t += det_t;
		//		if(count%10==0)
		std::cout << "t= " << t << ":" << element[N1 - 1].flux << "  " << element[N1].flowIn << " " << element[N1 + N2].flowIn << std::endl;
		updateElement(elementDown, t, 1);
		//		coutChangedData(element);
		//		if (count % 10 == 0)
		{
			pltTSDensity2FileforUnity(elementDown, timeInt);
			pltTSDensity2FileforTecplot(elementDown, timeInt, t);
		}
	}
	system("pause");
	return 0;
}