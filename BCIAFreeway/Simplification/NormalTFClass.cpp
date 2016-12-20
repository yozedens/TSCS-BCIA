#include"NormalTFClass.h"

Element::Element()
{
	density = flowIn = flux = omega = 0;
	EOut[0] = EOut[1] = -1;
	flag = 10;
	rampIn = rampOut = 0;
	laneN = 3;
	X = Y = 0;
}

Element::Element(double x, double y)
{
	X = x;
	Y = y;
}
Element::~Element()
{
}

