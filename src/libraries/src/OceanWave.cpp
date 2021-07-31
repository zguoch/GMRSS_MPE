

#include "proj.h"
#include "OceanWave.h"

// 服从均匀分布的随机数生成函数
double uniform(double a, double b)
{
	return (a+(b-a)*rand()/RAND_MAX);
}

// 服从正态分布的随机数生成函数
double Gaussrand(double Ave, double SD)
{
	static double U,V;
	static int phase=0;
	double Z;

	if(phase==0)
	{
		U=rand()/(RAND_MAX+1.0);
		V=rand()/(RAND_MAX+1.0);
		Z=sqrt(-2.0*log(U))*sin(2.0*PI*V);
	}
	else
	{
		Z=sqrt(-2.0*log(U))*cos(2.0*PI*V);
	}

	phase=1-phase;
	Z=Ave+Z*SD;
	return Z;
}

// P-M海浪谱函数
double PM1(double w,double U19p5,double gravity)
{
	return 8.1/1000.0*gravity*gravity/pow(w,5.0)*exp(-0.74*pow(gravity/(w*U19p5),4.0));
}

// 二维P-M海浪谱函数
double PM2(double w,double U19p5,double theta,double gravity)
{
	return 8.1/1000.0*gravity*gravity/pow(w,5.0)*exp(-0.74*pow(gravity/(w*U19p5),4.0))*8.0/(3.0*PI)*pow(cos(theta),4.0);
}

// 10米高度风速计算任意高度风速的函数
double Wind(double z,double U10)
{
	return U10*(1.0+sqrt(0.5*sqrt(U10)/1000.0)/0.4*log(z/10.0));
}