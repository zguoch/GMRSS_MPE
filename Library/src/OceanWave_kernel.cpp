// 为了方便在desktop和cmd中同时调用，二者的区别只有函数形参和进度条
#ifdef OCEANWAVE_CMD
    double WaveHeight(const Par_OceanWave& parm, std::vector<std::vector<std::vector<double> > >& h)
    {
#elif OCEANWAVE_DESKTOP_APP
    int MainWindow::doOceanWave()
    {
        OCEANWAVE::Par_OceanWave parm = m_par_OceanWave;
#endif
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
    double U19p5=OCEANWAVE::Wind(19.5,U10);
    // 生成0-2PI之间的随机数
    long i,j;
    double** num = new double*[M];
    num[0] = new double[M*N];
    for(i=1;i<M;++i)
    {
        num[i]=num[i-1]+N;		
    }
    if(parm.useSrand)srand(time(NULL));
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            num[i][j]=OCEANWAVE::uniform(0.0,2.0*PI);
        }
    }
    
    long ix,iy,iz,it,m,n;
    // long nProcess=0;
    double x,y,z,t,sum,sumx,sumy,sumz,Vx,Vy,Vz;
    // double k,omega,theta;

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
    // 先计算omega, k, theta并存储到变量里
    // double k,omega,theta;
    std::vector<double> arr_omega(M), arr_k(M), arr_theta(N);
    std::vector<std::vector<double> > factor(M);
    for (size_t i = 0; i < M; i++)factor[i].resize(N);
    
    for(m=0;m<M;m++)
    {
        arr_omega[m]=omega1+m*d_omega;
        arr_k[m]=arr_omega[m]*arr_omega[m]/gravity;
        for(n=0;n<N;n++)
        {
            arr_theta[n]=theta1+n*d_theta;	
            factor[m][n]=sqrt(2.0*OCEANWAVE::PM2(arr_omega[m],U19p5,arr_theta[n]-alfa,gravity)*d_omega*d_theta);
            // sum+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
        }	
    }
    omp_set_num_threads(parm.nThreads);
    #ifdef OCEANWAVE_CMD
        MultiProgressBar multibar_h(Ny*Nt,COLOR_BAR_BLUE);
    #elif OCEANWAVE_DESKTOP_APP
        ui->roundProgressBar->setRange(0, Ny*Nt);
        int ind=0;
    #endif
    for(it=0;it<Nt;it++)
    {
        double t=t1+it*d_t;
        #pragma omp parallel for private(iy, ix, m,n, y, x, sum) //omega, k, theta, 
        for(iy=0;iy<Ny;iy++)
        {
            y=y1+iy*d_y;
            for(ix=0;ix<Nx;ix++)
            {
                x=x1+ix*d_x;
                sum=0.0;
                
                for(m=0;m<M;m++)
                {
                    // omega=omega1+m*d_omega;
                    // k=omega*omega/gravity;
                    for(n=0;n<N;n++)
                    {
                        // theta=theta1+n*d_theta;	 
                        // sum+=sqrt(2.0*PM2(arr_omega[m],U19p5,arr_theta[n]-alfa,gravity)*d_omega*d_theta)*cos(arr_k[m]*x*cos(arr_theta[n])+arr_k[m]*y*sin(arr_theta[n])-arr_omega[m]*t+num[m][n]);
                        sum+=factor[m][n]*cos(arr_k[m]*x*cos(arr_theta[n])+arr_k[m]*y*sin(arr_theta[n])-arr_omega[m]*t+num[m][n]);
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
            #ifdef OCEANWAVE_CMD
                if(parm.showProgress){
                    #pragma omp critical
                    multibar_h.Update();
                }
            #elif OCEANWAVE_DESKTOP_APP
                if(parm.showProgress){
                    #pragma omp critical
                    ind++;
                    ui->roundProgressBar->setValue(ind);
                }
            #endif
            
        }
        // save result at t
        SaveResult(parm,wave_h, wave_U, t, SAVE_RESULT_h);
    }

    // 开始计算模拟海水速度数据
    #ifdef OCEANWAVE_CMD
        MultiProgressBar multibar_U(Nz*Nt*Ny,COLOR_BAR_BLUE);
    #elif OCEANWAVE_DESKTOP_APP
        ind = 0;
        ui->roundProgressBar->setRange(0, Nz*Nt*Ny);
        ui->roundProgressBar->setValue(1);
    #endif
    
    for(it=0;it<Nt;it++)
    {
        t=t1+it*d_t;
        for(iz=0;iz<Nz;iz++)
        {
            z=z1+iz*d_z;
            #pragma omp parallel for private(iy, ix, m, n, y, x, sumx, sumy, sumz, Vx, Vy, Vz) //omega, k, theta, 
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
                        // omega=omega1+m*d_omega;
                        // k=omega*omega/gravity;
                        for(n=0;n<N;n++)
                        {
                            // theta=theta1+n*d_theta;	
                            // sumx+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
                            // sumy+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*sin(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
                            // sumz+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
                            sumx+=factor[m][n]*exp(-arr_k[m]*z)*arr_omega[m]*cos(arr_theta[n])*sin(arr_k[m]*x*cos(arr_theta[n])+arr_k[m]*y*sin(arr_theta[n])-arr_omega[m]*t+num[m][n]);
                            sumy+=factor[m][n]*exp(-arr_k[m]*z)*arr_omega[m]*sin(arr_theta[n])*sin(arr_k[m]*x*cos(arr_theta[n])+arr_k[m]*y*sin(arr_theta[n])-arr_omega[m]*t+num[m][n]);
                            sumz+=factor[m][n]*exp(-arr_k[m]*z)*arr_omega[m]*cos(arr_k[m]*x*cos(arr_theta[n])+arr_k[m]*y*sin(arr_theta[n])-arr_omega[m]*t+num[m][n]);
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
                #ifdef OCEANWAVE_CMD
                    if(parm.showProgress){
                    #pragma omp critical
                    multibar_U.Update();
                    }
                #elif OCEANWAVE_DESKTOP_APP
                    if(parm.showProgress){
                        #pragma omp critical
                        ind++;
                        ui->roundProgressBar->setValue(ind);
                    }
                #endif
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
