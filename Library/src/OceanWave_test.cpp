/**
 * @file OceanWave_test.cpp
 * @author 杜劲松，郭志馗 (jinsongdu@cug.edu.cn, zguo@geomar.de)
 * @brief 测试二维海浪波高与海水流速模拟
 * @version 0.1
 * @date 2021-07-31
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "proj.h"
#include "OceanWave.h"
using namespace OCEANWAVE;
int main()
{
	time_t start,end;

    // 输入参数
	// 重力加速度(m/s^2)
	double gravity=9.81;

	// X方向空间范围和采样间隔(m)
	double x1=-50.0;
	double x2=+50.0;
	double d_x=1.0;

	// Y方向空间范围和采样间隔(m)
	double y1=-50.0;
	double y2=+50.0;
	double d_y=1.0;

	// Z方向空间范围和采样间隔(m)
	double z1=0.0;
	double z2=+10.0;
	double d_z=1.0;

	// 时间方向范围和采样间隔(s)
	double t1=0.0;
	double t2=+1000.0;
	double d_t=500.0;

	// 主浪方向(度)
	double alfa=45.0*PI/180.0;

	// 海面以上10米高度风速(m/s)
	double U10=20.0;

	// 最大角频率(rad/s)
	double omega2=3.0;

	// 角频率采样间隔(rad/s)
	double d_omega=0.01;

	// 方位采样间隔(度)
	double d_theta=1.0*PI/180.0;

	// 输入保存波高数据的文件名
	char DataOut1[256]="波高模拟数据.txt";

	// 输入保存海水速度数据的文件名
	char DataOut2[256]="海水速度模拟数据.txt";

    /////////////////////////////////////////////////////////////////////
	long Nx=(long)((x2-x1)/d_x+1.0);
	long Ny=(long)((y2-y1)/d_y+1.0);
	long Nz=(long)((z2-z1)/d_z+1.0);
	long Nt=(long)((t2-t1)/d_t+1.0);

	double omega1=d_omega;
	long M=(long)((omega2-omega1)/d_omega+1.0);

	double theta1=-90*PI/180.0;
	double theta2=+90.0*PI/180.0;
	long N=(long)((theta2-theta1)/d_theta+1.0);

	// 转换计算海面以上19.5米处的风速
	double U19p5=Wind(19.5,U10);

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
	double x,y,z,t,k,omega,sum,sumx,sumy,sumz,h,Vx,Vy,Vz,theta;

	// 打开保存波高数据的文件
	FILE *fpOut1=fopen(DataOut1, "w+");
	// 向波高数据文件输入文件头
	fprintf(fpOut1,"x[m]\ty[m]\tt[s]\th[m]\n");

	// 开始计算模拟波高数据
    printf("\nCalculating process : 1234");
	time(&start);
	for(it=0;it<Nt;it++)
    {
		nProcess=(long)((it+1)*100/Nt);
		printf("\b\b\b\b%2ld%% ",nProcess);

		t=t1+it*d_t;
		for(iy=0;iy<Ny;iy++)
		{
			y=y1+iy*d_y;
			for(ix=0;ix<Nx;ix++)
			{
				x=x1+ix*d_x;
				sum=0.0;
				for(m=0;m<M;m++)
				{
					omega=omega1+m*d_omega;
					k=omega*omega/gravity;
					for(n=0;n<N;n++)
					{
						theta=theta1+n*d_theta;	
						sum+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
					}	
				}
				h=sum;
				// 输出计算结果
				fprintf(fpOut1,"%.3lf\t%.3lf\t%.3lf\t%.3lf\n",x,y,t,h);
			}
		}
	}
	fclose(fpOut1);
	printf("\n\nThe wave height data has been simulated!\n\n");

	// 打开保存海水速度数据的文件
	FILE *fpOut2=fopen(DataOut2, "w+");
	// 向海水速度数据文件输入文件头
	fprintf(fpOut2,"x[m]\ty[m]\tz[m]\tt[s]\tVx[m/s]\tVy[m/s]\tVz[m/s]\n");

	// 开始计算模拟海水速度数据
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
					sumx=0.0;
					sumy=0.0;
					sumz=0.0;
					for(m=0;m<M;m++)
					{
						omega=omega1+m*d_omega;
						k=omega*omega/gravity;
						for(n=0;n<N;n++)
						{
							theta=theta1+n*d_theta;	
							sumx+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
							sumy+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*sin(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
							sumz+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
						}	
					}
					Vx=sumx;
					Vy=sumy;
					Vz=sumz;
					// 输出计算结果
					fprintf(fpOut2,"%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",x,y,z,t,Vx,Vy,Vz);
				}
			}
		}
	}
	fclose(fpOut2);
	printf("\n\nThe velocity data has been simulated!\n\n");

	delete[] num[0]; 
	delete[] num;
	
	time(&end);
	printf("\n\nCalculation is finished !");
	printf("\n\nComputing time: %.5lf seconds",difftime(end,start));
   
	// 计算结束
    printf("\n\nCalculating completed!\t\t\n\nPlease press any key to exit ...\n\n");
    // getch();
    return 0;
}
