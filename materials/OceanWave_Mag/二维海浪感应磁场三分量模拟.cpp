// 头文件
#include <stdio.h>
#include <math.h>
#include <mbstring.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <memory.h>

// 圆周率
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

// 服从均匀分布的随机数生成函数
double uniform(double a, double b);
// 服从正态分布的随机数生成函数
double Gaussrand(double Ave, double SD);
// P-M海浪谱函数(圆频率)
double PM1(double w,double U19p5,double gravity);
// P-M海浪谱函数(角波数)
double PM2(double k,double U19p5,double gravity);
// 二维P-M海浪谱&ISSc方向谱函数(圆频率)
double PM3(double w,double theta,double U19p5,double gravity);
// 波数转圆频率函数
double Convert_A(double k,double h,double gravity);
double Convert_B(double k,double gravity);
double Convert_C(double omega,double gravity);
// 10米高度风速计算任意高度风速的函数
double Wind(double h,double U10);

int main()
{
	time_t start,end;

    // 输入参数
	// 重力加速度(m/s^2)
	double gravity=9.81;

	// X方向空间范围和采样间隔(m)
	double x1=0.0;
	double x2=+0.0;
	double d_x=1.0;

	// Y方向空间范围和采样间隔(m)
	double y1=0.0;
	double y2=+0.0;
	double d_y=1.0;

	// Z方向空间范围和采样间隔(m)
	double z1=100.0;
	double z2=+100.0;
	double d_z=1.0;

	// 时间方向范围和采样间隔(s)
	double t1=-1000.0;
	double t2=+1000.0;
	double d_t=1.0;

	// 主浪方向(度)
	double alfa=45.0;

	// 海面以上10米高度风速(m/s)
	double U10=26.0;

	// 最大角频率(rad/s)
	double omega2=1.0;

	// 角频率采样间隔(rad/s)
	double d_omega=0.001;

	// 方位采样间隔(度)
	double d_theta=1.0;

	// 海水磁导率(H/m)
	double Miu=1.0;

	// 海水电导率(S/m)
	double Sigma=4.5;

	// 主磁场强度(nT)
	double F=44000.0;
	// 主磁场倾角(度)
	double I=20.0;
	// 主磁场偏角(度)
	double D=10.0;

	// 输入保存海浪磁场响应模拟数据的文件名
	char DataOut[256]="海浪磁场响应模拟数据.txt";

    /////////////////////////////////////////////////////////////////////
	long Nx=(long)((x2-x1)/d_x+1.0);
	long Ny=(long)((y2-y1)/d_y+1.0);
	long Nz=(long)((z2-z1)/d_z+1.0);
	long Nt=(long)((t2-t1)/d_t+1.0);

	double omega1=d_omega;
	long M=(long)((omega2-omega1)/d_omega+1.0);

	double theta1=-90*PI/180.0;
	double theta2=+90.0*PI/180.0;
	d_theta=d_theta*PI/180.0;
	long N=(long)((theta2-theta1)/d_theta+1.0);

	alfa=alfa*PI/180.0;

	Miu=Miu*4.0*PI*1E-007;

	// 转换计算海面以上19.5米处的风速
	double U19p5=Wind(19.5,U10);
	printf("\n%.5lf\n",U19p5);

	I=I*PI/180.0;
	D=D*PI/180.0;
	double Fx=F*cos(I)*cos(D);
	double Fy=F*cos(I)*sin(D);
	double Fz=F*sin(I);
	
	// 生成0-2PI之间的随机数
	long i,j;
	double** num = new double*[M];
	num[0] = new double[M*N];
	for(i=1;i<M;++i)
    {
        num[i]=num[i-1]+N;		
    }
	srand(time(NULL));
	for(i=0;i<M;i++)
	{
		for(j=0;j<N;j++)
		{
			num[i][j]=uniform(0.0,2.0*PI);
		}
	}
	
	long ix,iy,iz,it,m,n;
	long nProcess=0;
	double x,y,z,t,k,omega,sum,sum_Bx,sum_By,sum_Bz,Bx,By,Bz,theta,fai1,fai2;

	// 打开保存海浪磁场响应模拟数据的文件
	FILE *fpOut=fopen(DataOut, "w+");
	// 向波高数据文件输入文件头
	fprintf(fpOut,"x[m]\ty[m]\tz[m]\tt[s]\tBx[nT]\tBy[nT]\tBz[nT]\n");

	// 开始计算海浪磁场响应模拟数据
    printf("\nCalculating process : 1234");
	time(&start);
	for(it=0;it<Nt;it++)
    {
		nProcess=(long)((it+1)*100/Nt);
		printf("\b\b\b\b%2ld%% ",nProcess);

		t=t1+it*d_t;
		for(iz=0;iz<Nz;iz++)
		{
			z=z1+iz*d_z;
			for(iy=0;iy<Ny;iy++)
			{
				y=y1+iy*d_y;
				for(ix=0;ix<Nx;ix++)
				{
					x=x1+ix*d_x;
					sum_Bx=0.0;
					sum_By=0.0;
					sum_Bz=0.0;
					for(m=0;m<M;m++)
					{
						omega=omega1+m*d_omega;
						k=omega*omega/gravity;
						for(n=0;n<N;n++)
						{
							theta=theta1+n*d_theta;	
							fai1=atan2(Fx*cos(theta)+Fy*sin(theta),Fz);
							fai2=PI/2.0-fai1;
							if(z>=0.0)
							{
								sum=sqrt(2.0*PM3(omega,theta-alfa,U19p5,gravity)*d_omega*d_theta)*omega/k*exp(-k*z)*sqrt(Fx*Fx*cos(theta)*cos(theta)+Fy*Fy*sin(theta)*sin(theta)+Fz*Fz);
								sum_Bx+=sum*(-2.0*k*z+1.0)*cos(theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+fai1+num[m][n]);
								sum_By+=sum*(-2.0*k*z+1.0)*sin(theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+fai1+num[m][n]);
								sum_Bz+=sum*(+2.0*k*z+1.0)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t-fai2+num[m][n]);
							}
							else
							{
								sum=sqrt(2.0*PM3(omega,theta-alfa,U19p5,gravity)*d_omega*d_theta)*omega/k*exp(k*z)*sqrt(Fx*Fx*cos(theta)*cos(theta)+Fy*Fy*sin(theta)*sin(theta)+Fz*Fz);
								sum_Bx+=sum*cos(theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+fai1+num[m][n]);
								sum_By+=sum*sin(theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+fai1+num[m][n]);
								sum_Bz+=sum*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t-fai2+num[m][n]);
							}	
						}	
					}
					Bx=sum_Bx*Miu*Sigma/4.0;
					By=sum_By*Miu*Sigma/4.0;
					Bz=-sum_Bz*Miu*Sigma/4.0;
					// 输出计算结果
					fprintf(fpOut,"%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.15lf\t%.15lf\t%.15lf\n",x,y,z,t,Bx,By,Bz);
				}
			}
		}
	}
	delete[] num[0]; 
	delete[] num;
	fclose(fpOut);
	time(&end);
	printf("\n\nCalculation is finished !");
	printf("\n\nComputing time: %.5lf seconds",difftime(end,start));
   
	// 计算结束
    printf("\n\nCalculating completed!\t\t\n\nPlease press any key to exit ...\n\n");
    getch();
    return 0;
}

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

// P-M海浪谱函数(圆频率)
double PM1(double w,double U19p5,double gravity)
{
	return 8.1/1000.0*gravity*gravity/pow(w,5.0)*exp(-0.74*pow(gravity/(w*U19p5),4.0));
}
// P-M海浪谱函数(角波数)
double PM2(double k,double U19p5,double gravity)
{
	return 8.1/1000.0/pow(k,3.0)*exp(-0.74*pow(gravity,2.0)/(pow(k,2.0)*pow(U19p5,4.0)))/2.0;
}
// 二维P-M海浪谱&ISSc方向谱函数(圆频率)
double PM3(double w,double theta,double U19p5,double gravity)
{
	return 8.1/1000.0*gravity*gravity/pow(w,5.0)*exp(-0.74*pow(gravity/(w*U19p5),4.0))*8.0*pow(cos(theta),4.0)/(3.0*PI);
}

// 波数转圆频率函数
double Convert_A(double k, double h,double gravity)
{
	return sqrt(k*gravity*(1+k*k/(363.0*363.0))*tanh(k*h));

}
double Convert_B(double k,double gravity)
{
	return sqrt(k*gravity);
}
double Convert_C(double omega,double gravity)
{
	return omega*omega/gravity;
}

// 10米高度风速计算任意高度风速的函数
double Wind(double h,double U10)
{
	return U10*(1.0+sqrt(0.5*sqrt(U10)/1000.0)/0.4*log(h/10.0));
}