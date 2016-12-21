#pragma once
#include"NormalTFClass.h"
#ifndef SETVALUE_H
#define SETVALUE_H

const double PI = 3.14159265358979323846;//预定义圆周率常量
const double alpha1 = 50.0 / 180.0*PI, alpha2 = 25.0 / 180.0*PI, alpha3 = 55.0 / 180.0*PI, alpha4 = 95.0 / 180.0*PI;//与水平方向的逆时针到角

const double T = 400;//时间长度
const double Det_x = 20;//网格长度
extern double det_t;//可变时间间隔 

const double Length1 = 15600;//第一支线长度
const int N1 = 780;//第一支线网格数
const double Length11 = 900;//第一支线第一段长度
const int N11 = 45;//第一支线第一段网格数
const double Length12 = 10000;//第一支线第二段长度
const int N12 = 500;//第一支线第二段网格数
const double Length13 = 4700;//第一支线第三段长度
const int N13 = 235;//第一支线第三段网格数

const double Length2 = 5200;//第二支线长度
const int N2 = 260;//第二支线网格数
const double Length21 = 2800;//第二支线第一段长度
const int N21 = 140;//第二支线第一段网格数
const double Length22 = 2400;//第二支线第二段长度
const int N22 = 120;//第二支线第二段网格数

const int N = 1040;//总网格数
const int NLanes = 3;//单向车道数

const double vMax = 20;//最大速度
const double densityJam = 0.16;//最大密度
const double densityCr = densityJam / 2.0;//临界密度
const double flowCr = vMax*densityJam / 4.0;//最大流量
const double RateFlowDiverge = 0.4;//分岔口分流比率
const double boundryDensity = 0.2*densityJam;//进口处固定密度
//const double boundryVelocity = 0.4*vMax;//进口处固定速度

const int NI = 13;//上行入口数目
const int NO = 12;//上行出口数目
const double NIx[NI] = {    1740,     2000,2300,4020,4700,7100,  8600,11000,   13900,      16700,17000,18300,19100 };//上行入口位置
const double NOx[NO] = { 1400,   1900,2100,     4180,     6300,7600,  10800,13300,   15100,16700,17800,          20200};//上行出口位置

const int NID = 7;//下行入口数目
const int NOD = 11;//下行出口数目
const double NIDx[NID] = {                         6980,  10600,  13700,  15000, 16400,  17800,   20300};//下行入口位置
const double NODx[NOD] = { 1500,2100,2450,4800,6410,  8800,  11200,  14000,        16900,  19000,19300 };//下行出口位置

const double tau = 30;//延时
const double beta = tau*vMax / densityJam;//系数

double CFL_Det_t(Element*);//齐次方程CFL条件
double CFL_Det_t_ST(Element*);//非齐次方程CFL条件
double CFL_Det_t_yzd(Element*);//yzd自拟定CFL条件
double iniOmega(double);

void ini(Element*, int UorD = 0);//全局初始化函数,UorD=0表示上行，=1表示下行
void updateEIO(Element*, int UorD = 0);//上下游网格更新函数

void updateFlux(Element*);//更新数值流通量函数
double cal_g(double t);//自定义时间函数
void Fb(Element*, int UorD = 0);//特殊流量控制
void updateST(Element*, double t, int UorD = 0);//流入流出更改函数，模拟匝道io
void updateElement(Element*, double, int UorD = 0);//网格密度更新函数
double FGodunov(double, double);//Godunov数值流通量

void trans(double &x, double &y, double alpha);//向量逆时针旋转alpha弧度

double veLaneN(double, int);
double VLaneN(double, int);
double flowLaneN(double, int);
double FGLaneN(double, double, int);
double demandN(double, int);
double supplyN(double, int);

#endif
