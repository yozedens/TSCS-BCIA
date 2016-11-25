#include<iostream>
#include"NormalTFClass.h"
#include"SetValue.h"
#include"pltData.h"
#include<time.h>
using namespace std;

int main()
{
	Element element[N];
	ini(element);
	double t = 0, count = 1;
	
	int timeInt = time(0);
//	pltTSDensity2File(element, timeInt);
//	pltTSDensity2File2(element, 0);
	pltTSDensity2File3(element, 0);
//	pltGridPosition2File(element);
//	pltGridPosition2FileforUnity(element);
	updateEIO(element);
	
	while (t < T)
	{		
		updateElement(element);
		if(int(t)%50==0)
			std::cout <<"t= "<<t <<":"<<element[N1 - 1].flux << "  " << element[N1].flowIn << " " << element[N1 + N2].flowIn << std::endl;
//		coutChangedData(element);
		t += det_t;
		pltTSDensity2File(element,timeInt); 
//		pltTSDensity2File2(element, t);
		pltTSDensity2File3(element, t);
		count++;
	}

	cout <<count <<" lines!\n";
	system("pause");
	return 0;
}