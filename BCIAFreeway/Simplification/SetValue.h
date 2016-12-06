#pragma once
#include"NormalTFClass.h"
#ifndef SETVALUE_H
#define SETVALUE_H

const double PI = 3.14159265358979323846;//预定义圆周率常量
const double alpha1 = 50.0 / 180.0*PI, alpha2 = 25.0 / 180.0*PI, alpha3 = 55.0 / 180.0*PI, alpha4 = 95.0 / 180.0*PI;//与水平方向的逆时针到角
const double alpha5 = 355.0 / 180.0*PI, alpha6 = 325.0 / 180.0*PI, alpha7 = 65.0 / 180.0*PI;//与水平方向的逆时针到角

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

const double Length3 = 4800;//第三支线长度
const int N3 = 240;//第三支线网格数
const double Length31 = 2000;//第三支线第一段长度
const int N31 = 100;//第三支线第一段网格数
const double Length32 = 1400;//第三支线第二段长度
const int N32 = 70;//第三支线第二段网格数
const double Length33 = 1400;//第三支线第三段长度
const int N33 = 70;//第三支线第三段网格数

const int N = 1280;//总网格数
const int NLanes = 3;//单向车道数

const double vMax = 20;//最大速度
const double densityJam = 0.16;//最大密度
const double densityCr = densityJam / 2.0;//临界密度
const double flowCr = vMax*densityJam / 4.0;//最大流量
const double RateFlowDiverge = 0.4;//分岔口分流比率

const int NI = 12;//单边入口数目
const int NO = 10;//单边出口数目
const double NIx[NI] = { 1740,     2150,2300,4150,4700,7100,  8600,11000,   13900,18650,19100,24250 };//单边入口位置/*大山桥处前移400*//*北皋桥处前移100*//*华谊桥处前移250*/
const double NOx[NO] = { 1400,1900,2100,     4100,     6300,7600,  10800,13300,   18200,      24200 };//单边出口位置/*大山桥出后移400*/0/*北皋桥出后移100*//*华谊桥处后移200*/

const double tau = 30;//延时
const double beta = tau*vMax / densityJam;//系数

double ve(double);//速密关系函数
double V(double);
double flow(double);//流量密度函数
double demand(double);//需求函数
double supply(double);//供给函数
double CFL_Det_t(Element*);//齐次方程CFL条件
double CFL_Det_t_ST(Element*);//非齐次方程CFL条件
double CFL_Det_t_yzd(Element*);//yzd自拟定CFL条件
double intOmega(double);

void ini(Element*);//全局初始化函数
void updateFlux(Element*);//更新数值流通量函数
double cal_g(double t);//自定义时间函数
void updateST(Element*, double t);//流入流出更改函数，模拟匝道io
void updateEIO(Element*);//上下游网格更新函数
void updateElement(Element*, double);//网格密度更新函数
double FGodunov(double, double);//Godunov数值流通量

void trans(double &x, double &y, double alpha);//向量逆时针旋转alpha弧度

double veLaneN(double, int);
double VLaneN(double, int);
double flowLaneN(double, int);
double FGLaneN(double, double, int);
double demandN(double, int);
double supplyN(double, int);

#endif
