#pragma once
class Element
{
public:
	Element();
	Element(double, double);
	double density, velocity, flowIn, flux;
	double omega;
	double rampIn, rampOut;
	double laneN;
	~Element();
	void initialization(double a)
	{
		density = a;
	}
	int EOut[2];
	int flag;
	void setValueXY(double x, double y) { X = x; Y = y; }
	double getValueX() { return X; }
	double getValueY() { return Y; }
private:
	double X, Y;

};
