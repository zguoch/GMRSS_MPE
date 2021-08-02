 

#include "proj.h"
#include "OceanWave.h"
#include "MultiProgressBar.h"
#include "omp.h"
#include <fstream>
#include <cstdio>
namespace OCEANWAVE
{
	// void Par_OceanWave::print()
	

	double WaveHeight_origin(const Par_OceanWave& parm, std::vector<std::vector<std::vector<double> > >& h)
	{
		time_t start,end;
		time(&start);
		// 开始计算之前，先判断输出文件.txt格式是否已经存在，如果存在则删除，以为它是以app的形式打开的，会续写在老的文件里面
		{
			std::string fname_h = parm.fname_WaveHeight+"."+parm.fmt_outputFile;
			std::string fname_U = parm.fname_SeawaterVelocity+"."+parm.fmt_outputFile;
			std::remove(fname_h.c_str());
			std::remove(fname_U.c_str());
		}
		// 参数
		double gravity=parm.gravity;
		// X方向空间范围和采样间隔(m)
		double x1=parm.minDmax_xyzt[0][0];
		double x2=parm.minDmax_xyzt[0][2];
		double d_x=parm.minDmax_xyzt[0][1];
		// Y方向空间范围和采样间隔(m)
		double y1=parm.minDmax_xyzt[1][0];
		double y2=parm.minDmax_xyzt[1][2];
		double d_y=parm.minDmax_xyzt[1][1];
		// Z方向空间范围和采样间隔(m)
		double z1=parm.minDmax_xyzt[2][0];
		double z2=parm.minDmax_xyzt[2][2];
		double d_z=parm.minDmax_xyzt[2][1];
		// 时间方向范围和采样间隔(s)
		double t1=parm.minDmax_xyzt[3][0];
		double t2=parm.minDmax_xyzt[3][2];
		double d_t=parm.minDmax_xyzt[3][1];
		// 主浪方向(度)
		double alfa=parm.alpha*PI/180.0;
		// 海面以上10米高度风速(m/s)
		double U10=parm.U10;
		// 最大角频率(rad/s)
		double omega2=parm.omega2;
		// 角频率采样间隔(rad/s)
		double d_omega=parm.dOmega;
		// 方位采样间隔(度)
		double d_theta=parm.dTheta*PI/180.0;

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
		if(parm.useSrand)srand(time(NULL)); //用户可选择是否使用srand
		for(i=0;i<M;i++)
		{
			for(j=0;j<N;j++)
			{
				num[i][j]=uniform(0.0,2.0*PI);
			}
		}
		
		long ix,iy,iz,it,m,n;
		// long nProcess=0;
		double x,y,z,t,k,omega,sum,sumx,sumy,sumz,Vx,Vy,Vz,theta;

		// // 根据传入的vector引用的大小判断是否返回数据，为什么要这么做？答：有时候如果用户给定的模拟数据量太大，就不用返回所有数据了，只写入文件就行
		// // 有这三种情况：（1）返回所有数据，也就是这个h是一个NtxNyxNx的矩阵；（2）只返回最后时刻的计算结果，h是一个1xNyxNx;(3)不返还任何数据，h的大小为0
		// int returnResults = RETURN_RESULTS_NO;
		// if(h.size()==0){
		// 	returnResults = RETURN_RESULTS_NO;
		// }else{
		// 	if(h.size()==1 && h[0].size()==Ny && h[0][0].size()==Nx){
		// 		returnResults = RETURN_RESULTS_LATEST;
		// 	}else if(h.size()==Nt && h[0].size()==Ny && h[0][0].size()==Nx){
		// 		returnResults = RETURN_RESULTS_ALL;
		// 	}else{
		// 		returnResults = RETURN_RESULTS_NO;
		// 	}
		// }
		// 对每一个时刻的结果进行保存
		std::vector<std::vector<double> > wave_h(Ny);
		for (size_t i = 0; i < Ny; i++)wave_h[i].resize(Nx);
		std::vector<std::vector<std::vector<std::vector<double> > > > wave_U(3);
		for (size_t i = 0; i < 3; i++){
			wave_U[i].resize(Nz);
			for (size_t iz = 0; iz < Nz; iz++){
				wave_U[i][iz].resize(Ny);
				for (size_t iy = 0; iy < Ny; iy++){
					wave_U[i][iz][iy].resize(Nx);
				}
			}
		}
		
		// 开始计算模拟波高数据
		omp_set_num_threads(parm.nThreads);
		MultiProgressBar multibar_U(Ny*Nt,COLOR_BAR_BLUE);
		for(it=0;it<Nt;it++)
		{
			double t=t1+it*d_t;
			#pragma omp parallel for private(iy, ix, m,n, omega, k, theta, y, x, sum)
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
					// // 根据情况判断是否需要把数据返回
					// if(returnResults == RETURN_RESULTS_LATEST && it==(Nt-1))
					// {
					// 	h[0][iy][ix]=sum;
					// }else if(returnResults == RETURN_RESULTS_ALL)
					// {
					// 	h[it][iy][ix]=sum;
					// }
					wave_h[iy][ix] = sum;
				}
				if(parm.showProgress){
					#pragma omp critical
					multibar_U.Update();
				}
			}
			// save result at t
			SaveResult(parm,wave_h, wave_U, t, SAVE_RESULT_h);
		}

		// 开始计算模拟海水速度数据
		MultiProgressBar multibar(Nz*Nt*Ny,COLOR_BAR_BLUE);
		for(it=0;it<Nt;it++)
		{
			t=t1+it*d_t;
			for(iz=0;iz<Nz;iz++)
			{
				z=z1+iz*d_z;
				#pragma omp parallel for private(iy, ix, m, n, omega, k, theta, y, x, sumx, sumy, sumz, Vx, Vy, Vz)
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
						// fprintf(fpOut2,"%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",x,y,z,t,Vx,Vy,Vz);
						wave_U[0][iz][iy][ix] = Vx;
						wave_U[1][iz][iy][ix] = Vy;
						wave_U[2][iz][iy][ix] = Vz;
					}
					if(parm.showProgress){
					#pragma omp critical
					multibar.Update();
				}
				}
			}
			SaveResult(parm,wave_h, wave_U, t, SAVE_RESULT_U);
		}
		//release pointer
		delete[] num[0]; 
		delete[] num;
		
		time(&end);
		printf("\n计算完毕 !");
		printf("\n用时: %.5lf 秒\n",difftime(end,start));

		return 0;
	}
	#define OCEANWAVE_CMD 1
	#include "OceanWave_kernel.cpp"
	void SaveResult(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h, 
		const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U)
	{
		if (parm.fmt_outputFile == "txt")
		{
			SaveResult_txt(parm, h, U, t, h_or_U);
		}else if(parm.fmt_outputFile == "grd")
		{
			SaveResult_grd(parm, h, U, t, h_or_U);
		}else if(parm.fmt_outputFile == "vtk")
		{
			SaveResult_VTK(parm, h, U, t, h_or_U);
		}
		else
		{
			std::cout<<"Error: 输出文件不支持"<<parm.fmt_outputFile<<"格式，请通过-f 指定文件格式，比如-f txt\n";
			exit(0);
		}
		
	}
	void SaveResult_VTK(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h, 
		const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U)
	{
		char str_t[256];
		sprintf(str_t, "%.0f", t);
		std::string fname_waveheight = (h_or_U==SAVE_RESULT_h ? parm.fname_WaveHeight : parm.fname_SeawaterVelocity)+"_"+str_t+"."+parm.fmt_outputFile;
		std::ofstream fout_waveheight(fname_waveheight);
		if(!fout_waveheight)
		{
			std::cout<<"Error: 打开文件失败, "<<fname_waveheight<<std::endl;
			exit(0);
		}
		fout_waveheight<<"# vtk DataFile Version 2.0\n";
		fout_waveheight<<"Result data generated by G&MRSSofMPE\n";
		fout_waveheight<<"ASCII\n";
		fout_waveheight<<"DATASET RECTILINEAR_GRID\n";
		int nx = (h_or_U==SAVE_RESULT_h ? h[0].size() : U[0][0][0].size());
		int ny = (h_or_U==SAVE_RESULT_h ? h.size() : U[0][0].size());
		int nz = (h_or_U==SAVE_RESULT_h ? 1 : U[0].size());
		fout_waveheight<<"DIMENSIONS "<<nx<<" "<<ny<<" "<<nz<<"\n";
		fout_waveheight<<"X_COORDINATES "<<nx<<" float\n";
		for (size_t ix = 0; ix < nx; ix++)fout_waveheight<<parm.minDmax_xyzt[0][0]+ix*parm.minDmax_xyzt[0][1]<<" ";
		fout_waveheight<<"\n";
		fout_waveheight<<"Y_COORDINATES "<<ny<<" float\n";
		for (size_t iy = 0; iy < ny; iy++)fout_waveheight<<parm.minDmax_xyzt[1][0]+iy*parm.minDmax_xyzt[1][1]<<" ";
		fout_waveheight<<"\n";
		fout_waveheight<<"Z_COORDINATES "<<nz<<" float\n";
		for (size_t iz = 0; iz < nz; iz++)fout_waveheight<<parm.minDmax_xyzt[2][0]+iz*parm.minDmax_xyzt[2][1]<<" ";
		fout_waveheight<<"\n";
		// wave height
		fout_waveheight<<"POINT_DATA "<<nx*ny*nz<<"\n";
		if(h_or_U == SAVE_RESULT_h)
		{
			fout_waveheight<<"SCALARS WaveHeight float\n";
			fout_waveheight<<"LOOKUP_TABLE default\n";
			for (size_t iy = 0; iy < ny; iy++)
			{
				for (size_t ix = 0; ix < nx; ix++)
				{
					fout_waveheight<<h[iy][ix]<<" ";
				}
				fout_waveheight<<"\n";
			}
		}else if(h_or_U == SAVE_RESULT_U)
		{
			fout_waveheight<<"VECTORS SeawaterVelocity float\n";
			// fout_waveheight<<"LOOKUP_TABLE default\n";
			for (size_t iz = 0; iz < nz; iz++)
			{
				for (size_t iy = 0; iy < ny; iy++)
				{
					for (size_t ix = 0; ix < nx; ix++)
					{
						fout_waveheight<<U[0][iz][iy][ix]<<" ";
						fout_waveheight<<U[1][iz][iy][ix]<<" ";
						fout_waveheight<<U[2][iz][iy][ix]<<"\n";
					}
					// fout_waveheight<<"\n";
				}
			}
		}
		fout_waveheight.close();
	}
	void SaveResult_grd(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h, 
		const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U)
	{
		// 速度不支持grd格式
		if(h_or_U == SAVE_RESULT_U)return;
		std::string fname_waveheight = (h_or_U==SAVE_RESULT_h ? parm.fname_WaveHeight : parm.fname_SeawaterVelocity)+"_"+std::to_string(t)+"."+parm.fmt_outputFile;
		std::ofstream fout_waveheight(fname_waveheight);
		if(!fout_waveheight)
		{
			std::cout<<"Error: 打开文件失败, "<<fname_waveheight<<std::endl;
			exit(0);
		}
		fout_waveheight<<"DSAA\n";
		fout_waveheight<<h[0].size()<<" "<<h.size()<<"\n";
		fout_waveheight<<parm.minDmax_xyzt[0][0]<<" "<<parm.minDmax_xyzt[0][2]<<"\n";
		fout_waveheight<<parm.minDmax_xyzt[1][0]<<" "<<parm.minDmax_xyzt[1][2]<<"\n";
		fout_waveheight<<0<<" "<<1<<"\n";
		for (size_t iy = 0; iy < h.size(); iy++)
		{
			for (size_t ix = 0; ix < h[0].size(); ix++)
			{
				fout_waveheight<<h[iy][ix]<<" ";
			}
			fout_waveheight<<"\n";
		}
		fout_waveheight.close();
	}
	void SaveResult_txt(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h, 
		const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U)
	{
		// .txt的情况下下把所有时间的数据写到一个文件里面，方便对单点随时间变化或剖面输出
		std::string filename = (h_or_U==SAVE_RESULT_h ? parm.fname_WaveHeight : parm.fname_SeawaterVelocity)+"."+parm.fmt_outputFile;
		std::ofstream fout_waveheight(filename, std::ofstream::app);
		if(!fout_waveheight)
		{
			std::cout<<"Error: 打开文件失败, "<<filename<<std::endl;
			exit(0);
		}
		if(h_or_U == SAVE_RESULT_h)
		{
			// 如果时间是最小之间，即第一个时刻则写入文件头
			if(t==parm.minDmax_xyzt[3][0])fout_waveheight<<"x[m]\ty[m]\tt[s]\th[m]\n";
			for (size_t iy = 0; iy < h.size(); iy++)
			{
				for (size_t ix = 0; ix < h[0].size(); ix++)
				{
					fout_waveheight<<parm.minDmax_xyzt[0][0]+ix*parm.minDmax_xyzt[0][1]<<"\t"
								<<parm.minDmax_xyzt[1][0]+iy*parm.minDmax_xyzt[1][1]<<"\t"
								<<t<<"\t"
								<<h[iy][ix]
								<<std::endl;
				}
				
			}
		}
		else if(h_or_U == SAVE_RESULT_U)
		{
			// 如果时间是最小之间，即第一个时刻则写入文件头
			if(t==parm.minDmax_xyzt[3][0])fout_waveheight<<"x[m]\ty[m]\tz[m]\tt[s]\tVx[m/s]\tVy[m/s]\tVz[m/s]\n";
			for (size_t iz = 0; iz < U[0].size(); iz++)
			{
				for (size_t iy = 0; iy < U[0][0].size(); iy++)
				{
					for (size_t ix = 0; ix < U[0][0][0].size(); ix++)
					{
						fout_waveheight<<parm.minDmax_xyzt[0][0]+ix*parm.minDmax_xyzt[0][1]<<"\t"
									   <<parm.minDmax_xyzt[1][0]+iy*parm.minDmax_xyzt[1][1]<<"\t"
									   <<parm.minDmax_xyzt[2][0]+iy*parm.minDmax_xyzt[2][1]<<"\t"
									   <<t<<"\t"
									   <<U[0][iz][iy][ix]<<"\t"
									   <<U[1][iz][iy][ix]<<"\t"
									   <<U[2][iz][iy][ix]<<"\t"
									   <<std::endl;
					}
				}
			}
		}
		fout_waveheight.close();
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
}
