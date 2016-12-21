#pragma once
#include"NormalTFClass.h"
#ifndef SETVALUE_H
#define SETVALUE_H

const double PI = 3.14159265358979323846;//Ԥ����Բ���ʳ���
const double alpha1 = 50.0 / 180.0*PI, alpha2 = 25.0 / 180.0*PI, alpha3 = 55.0 / 180.0*PI, alpha4 = 95.0 / 180.0*PI;//��ˮƽ�������ʱ�뵽��

const double T = 400;//ʱ�䳤��
const double Det_x = 20;//���񳤶�
extern double det_t;//�ɱ�ʱ���� 

const double Length1 = 15600;//��һ֧�߳���
const int N1 = 780;//��һ֧��������
const double Length11 = 900;//��һ֧�ߵ�һ�γ���
const int N11 = 45;//��һ֧�ߵ�һ��������
const double Length12 = 10000;//��һ֧�ߵڶ��γ���
const int N12 = 500;//��һ֧�ߵڶ���������
const double Length13 = 4700;//��һ֧�ߵ����γ���
const int N13 = 235;//��һ֧�ߵ�����������

const double Length2 = 5200;//�ڶ�֧�߳���
const int N2 = 260;//�ڶ�֧��������
const double Length21 = 2800;//�ڶ�֧�ߵ�һ�γ���
const int N21 = 140;//�ڶ�֧�ߵ�һ��������
const double Length22 = 2400;//�ڶ�֧�ߵڶ��γ���
const int N22 = 120;//�ڶ�֧�ߵڶ���������

const int N = 1040;//��������
const int NLanes = 3;//���򳵵���

const double vMax = 20;//����ٶ�
const double densityJam = 0.16;//����ܶ�
const double densityCr = densityJam / 2.0;//�ٽ��ܶ�
const double flowCr = vMax*densityJam / 4.0;//�������
const double RateFlowDiverge = 0.4;//�ֲ�ڷ�������
const double boundryDensity = 0.2*densityJam;//���ڴ��̶��ܶ�
//const double boundryVelocity = 0.4*vMax;//���ڴ��̶��ٶ�

const int NI = 13;//���������Ŀ
const int NO = 12;//���г�����Ŀ
const double NIx[NI] = {    1740,     2000,2300,4020,4700,7100,  8600,11000,   13900,      16700,17000,18300,19100 };//�������λ��
const double NOx[NO] = { 1400,   1900,2100,     4180,     6300,7600,  10800,13300,   15100,16700,17800,          20200};//���г���λ��

const int NID = 7;//���������Ŀ
const int NOD = 11;//���г�����Ŀ
const double NIDx[NID] = {                         6980,  10600,  13700,  15000, 16400,  17800,   20300};//�������λ��
const double NODx[NOD] = { 1500,2100,2450,4800,6410,  8800,  11200,  14000,        16900,  19000,19300 };//���г���λ��

const double tau = 30;//��ʱ
const double beta = tau*vMax / densityJam;//ϵ��

double CFL_Det_t(Element*);//��η���CFL����
double CFL_Det_t_ST(Element*);//����η���CFL����
double CFL_Det_t_yzd(Element*);//yzd���ⶨCFL����
double iniOmega(double);

void ini(Element*, int UorD = 0);//ȫ�ֳ�ʼ������,UorD=0��ʾ���У�=1��ʾ����
void updateEIO(Element*, int UorD = 0);//������������º���

void updateFlux(Element*);//������ֵ��ͨ������
double cal_g(double t);//�Զ���ʱ�亯��
void Fb(Element*, int UorD = 0);//������������
void updateST(Element*, double t, int UorD = 0);//�����������ĺ�����ģ���ѵ�io
void updateElement(Element*, double, int UorD = 0);//�����ܶȸ��º���
double FGodunov(double, double);//Godunov��ֵ��ͨ��

void trans(double &x, double &y, double alpha);//������ʱ����תalpha����

double veLaneN(double, int);
double VLaneN(double, int);
double flowLaneN(double, int);
double FGLaneN(double, double, int);
double demandN(double, int);
double supplyN(double, int);

#endif
