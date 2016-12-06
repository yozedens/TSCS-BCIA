#pragma once
#include"NormalTFClass.h"
#ifndef SETVALUE_H
#define SETVALUE_H

const double PI = 3.14159265358979323846;//Ԥ����Բ���ʳ���
const double alpha1 = 50.0 / 180.0*PI, alpha2 = 25.0 / 180.0*PI, alpha3 = 55.0 / 180.0*PI, alpha4 = 95.0 / 180.0*PI;//��ˮƽ�������ʱ�뵽��
const double alpha5 = 355.0 / 180.0*PI, alpha6 = 325.0 / 180.0*PI, alpha7 = 65.0 / 180.0*PI;//��ˮƽ�������ʱ�뵽��

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

const double Length3 = 4800;//����֧�߳���
const int N3 = 240;//����֧��������
const double Length31 = 2000;//����֧�ߵ�һ�γ���
const int N31 = 100;//����֧�ߵ�һ��������
const double Length32 = 1400;//����֧�ߵڶ��γ���
const int N32 = 70;//����֧�ߵڶ���������
const double Length33 = 1400;//����֧�ߵ����γ���
const int N33 = 70;//����֧�ߵ�����������

const int N = 1280;//��������
const int NLanes = 3;//���򳵵���

const double vMax = 20;//����ٶ�
const double densityJam = 0.16;//����ܶ�
const double densityCr = densityJam / 2.0;//�ٽ��ܶ�
const double flowCr = vMax*densityJam / 4.0;//�������
const double RateFlowDiverge = 0.4;//�ֲ�ڷ�������

const int NI = 12;//���������Ŀ
const int NO = 10;//���߳�����Ŀ
const double NIx[NI] = { 1740,     2150,2300,4150,4700,7100,  8600,11000,   13900,18650,19100,24250 };//�������λ��/*��ɽ�Ŵ�ǰ��400*//*�����Ŵ�ǰ��100*//*�����Ŵ�ǰ��250*/
const double NOx[NO] = { 1400,1900,2100,     4100,     6300,7600,  10800,13300,   18200,      24200 };//���߳���λ��/*��ɽ�ų�����400*/0/*�����ų�����100*//*�����Ŵ�����200*/

const double tau = 30;//��ʱ
const double beta = tau*vMax / densityJam;//ϵ��

double ve(double);//���ܹ�ϵ����
double V(double);
double flow(double);//�����ܶȺ���
double demand(double);//������
double supply(double);//��������
double CFL_Det_t(Element*);//��η���CFL����
double CFL_Det_t_ST(Element*);//����η���CFL����
double CFL_Det_t_yzd(Element*);//yzd���ⶨCFL����
double intOmega(double);

void ini(Element*);//ȫ�ֳ�ʼ������
void updateFlux(Element*);//������ֵ��ͨ������
double cal_g(double t);//�Զ���ʱ�亯��
void updateST(Element*, double t);//�����������ĺ�����ģ���ѵ�io
void updateEIO(Element*);//������������º���
void updateElement(Element*, double);//�����ܶȸ��º���
double FGodunov(double, double);//Godunov��ֵ��ͨ��

void trans(double &x, double &y, double alpha);//������ʱ����תalpha����

double veLaneN(double, int);
double VLaneN(double, int);
double flowLaneN(double, int);
double FGLaneN(double, double, int);
double demandN(double, int);
double supplyN(double, int);

#endif
