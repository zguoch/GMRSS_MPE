#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "GMRSS_MPEVersion.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    ,m_threadNumOMP(omp_get_max_threads())
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    statusBar()->showMessage(tr("G&MRSSofMPE version: ")+GMRSS_MPE_VERSION+"， 并行计算："+std::to_string(m_threadNumOMP).c_str()+"/"+std::to_string(omp_get_max_threads()).c_str()+" CPU");
    updateUILayout();
    initParameters();
    // thread for busy calculation
    watcher_ = new QFutureWatcher<int>;
    connect(watcher_, &QFutureWatcher<int>::finished,this, &MainWindow::busy_job_finished);

    // 暂时用不到右边的文本框，隐藏掉
    ui->splitter->setSizes(QList<int>() << 1 << 0);
    //round progress bar
    ui->roundProgressBar->setVisible(false);
    ui->roundProgressBar->setDecimals(0);
    QPalette pal = palette();
    // set background
    pal.setColor(QPalette::Background, Qt::white);
    pal.setColor(QPalette::Text, Qt::blue);
    pal.setBrush(QPalette::AlternateBase, Qt::white);
    ui->roundProgressBar->setAutoFillBackground(true);
    ui->roundProgressBar->setPalette(pal);
    QGradientStops gradientPoints;
    gradientPoints << QGradientStop(0, Qt::green) << QGradientStop(0.5, Qt::yellow) << QGradientStop(1, Qt::red);
    // and set it
    ui->roundProgressBar->setDataColors(gradientPoints);
    ui->roundProgressBar->setRange(0,100);
    ui->roundProgressBar->setValue(0);
    ui->roundProgressBar->setNullPosition(QRoundProgressBar::PositionBottom);

    ui->toolBar->setFixedHeight(36);
    ui->toolBar->setIconSize(QSize(36, 36));
    
    // ------------renderwindow------------
    m_renderWindow=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    this->ui->qvtkWidget2->SetRenderWindow(m_renderWindow);
    initRenderWindow();
    // ------2D chart---------
    // Set up the view
    m_vtkChartView=vtkSmartPointer<vtkContextView>::New();
    m_vtkChartView->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
    vtkSmartPointer<vtkChartXY> chart = 
        vtkSmartPointer<vtkChartXY>::New();
    m_vtkChartView->GetScene()->AddItem(chart);
    // 2. key for 2D charts
    m_vtkChartView->SetRenderWindow(this->ui->qvtkWidget->GetRenderWindow()); //must be here

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initRenderWindow()
{
    // Geometry
    vtkNew<vtkVectorText> text;
    // text->SetText("Gravity and Magnetic Response Simulation System of\nMarine Physical Environment\n(G&MRSSofMPE)\nVersion 1.0\n\nChina University of Geosciences (Wuhan)");
    text->SetText("G&MRSSofMPE V1.0");
    vtkNew<vtkElevationFilter> elevation;
    elevation->SetInputConnection(text->GetOutputPort());
    elevation->SetLowPoint(0,0,0);
    elevation->SetHighPoint(10,0,0);

    // Mapper
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(elevation->GetOutputPort());

    // Actor in scene
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);

    // VTK Renderer
    vtkNew<vtkRenderer> ren;

    // Add Actor to renderer
    ren->AddActor(actor);
    // ren->AddActor(test());

    ren->SetBackground2( 0.8, 0.4, 0.1 );
    ren->SetBackground( 0.1, 0.4, 0.8 );
    ren->TexturedBackgroundOn();
    this->ui->qvtkWidget2->GetRenderWindow()->AddRenderer(ren);
    this->ui->qvtkWidget2->GetRenderWindow()->Render();
}
void MainWindow::busy_job()
{
    ui->pushButton->setVisible(false);
    ui->roundProgressBar->setVisible(true);
    int index_tab=ui->tabWidget->currentIndex();
    switch (index_tab)
    {
    case 0:  //海浪的重磁响应
        {
            switch (getIndex_OceanWave_Grav_Mag())
            {
            case INDEX_OCEANWAVE_WAVE:
                {
                    getPar_OceanWave(m_par_OceanWave);
                    auto future = QtConcurrent::run(this, &MainWindow::doOceanWave);
                    watcher_->setFuture(future);
                }
                break;
            default:
                break;
            }
        }
        break;
    default:
        QMessageBox msgBox;
        msgBox.setWindowTitle(tr("Information"));
        msgBox.setText(ui->tabWidget->tabText(index_tab)+tr(": 正在开发中..."));
        msgBox.setStandardButtons(QMessageBox::Yes);
        msgBox.setDefaultButton(QMessageBox::Yes);
        msgBox.exec();
        //ui->textEdit->append("3D is comming soon");
        break;
    }
    
}
void MainWindow::Info(std::string infoText)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle(tr("Information"));
    msgBox.setText(tr(infoText.c_str()));
    msgBox.setStandardButtons(QMessageBox::Yes);
    msgBox.setDefaultButton(QMessageBox::Yes);
    msgBox.exec();
}
void MainWindow::busy_job_finished()
{
    // switch (m_dimension) {
    //     case 1:
    //     {
    //         //display
    //         ShowProps_1D();
    //     }
    //     break;
    //     case 2:
    //     {
    //         ShowProps_2D(ui->comboBox_selectProps->currentIndex(), m_xlabel, m_ylabel, m_zlabel,m_actorScale);
    //         m_vtkCameraInitialized=true;
    //     }
    //     break;
    // case 3:
    //     ui->textEdit->append(tr("3D is comming soon"));
    //     break;

    // }
    ui->roundProgressBar->setVisible(false); //calculation finished, hide progressbar
    ui->pushButton->setVisible(true);
}

void MainWindow::on_pushButton_clicked()
{
    busy_job();//using another thread to calculate, update chart/renderer after busy_job finishing
}
int MainWindow::getIndex_OceanWave_Grav_Mag()
{
    if(ui->radioButton_6->isChecked()){return INDEX_OCEANWAVE_WAVE;}
    else if(ui->radioButton_7->isChecked()) {return INDEX_OCEANWAVE_GRAV;}

    return -1;
}
void MainWindow::getPar_OceanWave(Par_OceanWave& par)
{
    par.xmin=ui->doubleSpinBox_xmin->value();
    par.dx=ui->doubleSpinBox_dx->value();
    par.xmax=ui->doubleSpinBox_xmax->value();
    par.ymin=ui->doubleSpinBox_ymin->value();
    par.dy=ui->doubleSpinBox_dy->value();
    par.ymax=ui->doubleSpinBox_ymax->value();
    par.zmin=ui->doubleSpinBox_zmin->value();
    par.dz=ui->doubleSpinBox_dz->value();
    par.zmax=ui->doubleSpinBox_zmax->value();
    
    par.alpha=ui->doubleSpinBox_alpha->value();
    par.dTheta=ui->doubleSpinBox_dTheta->value();
    par.U10=ui->doubleSpinBox_U10->value();
    par.omega2=ui->doubleSpinBox_omega2->value();
    par.dOmega=ui->doubleSpinBox_dOmega->value();
    par.tmin=ui->doubleSpinBox_tmin->value();
    par.dt=ui->doubleSpinBox_dt->value();
    par.tmax=ui->doubleSpinBox_tmax->value();
    par.gravity=ui->doubleSpinBox_gravity->value();

    par.fname_WaveHeight=ui->lineEdit->text().toStdString();
    par.fname_SeawaterVelocity=ui->lineEdit_2->text().toStdString();
}
int MainWindow::doOceanWave()
{
    time_t start,end;
    // 输入参数
	// 重力加速度(m/s^2)
	double gravity=m_par_OceanWave.gravity;
	// X方向空间范围和采样间隔(m)
	double x1=m_par_OceanWave.xmin;
	double x2=m_par_OceanWave.xmax;
	double d_x=m_par_OceanWave.dx;
	// Y方向空间范围和采样间隔(m)
	double y1=m_par_OceanWave.ymin;
	double y2=m_par_OceanWave.ymax;
	double d_y=m_par_OceanWave.dy;
	// Z方向空间范围和采样间隔(m)
	double z1=m_par_OceanWave.zmin;
	double z2=m_par_OceanWave.zmax;
	double d_z=m_par_OceanWave.dz;
	// 时间方向范围和采样间隔(s)
	double t1=m_par_OceanWave.tmin;
	double t2=m_par_OceanWave.tmax;
	double d_t=m_par_OceanWave.dt;
	// 主浪方向(度)
	double alfa=m_par_OceanWave.alpha*PI/180.0;
	// 海面以上10米高度风速(m/s)
	double U10=m_par_OceanWave.U10;
	// 最大角频率(rad/s)
	double omega2=m_par_OceanWave.omega2;
	// 角频率采样间隔(rad/s)
	double d_omega=m_par_OceanWave.dOmega;
	// 方位采样间隔(度)
	double d_theta=m_par_OceanWave.dTheta*PI/180.0;
	// // 输入保存波高数据的文件名
	// char DataOut1[256]=m_par_OceanWave.fname_WaveHeight.c_str();
	// // 输入保存海水速度数据的文件名
	// char DataOut2[256]=m_par_OceanWave.fname_SeawaterVelocity.c_str();

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
    std::ofstream fout1(m_par_OceanWave.fname_WaveHeight);
    if(!fout1)
    {
        std::cout<<"打开文件失败: "<<m_par_OceanWave.fname_WaveHeight<<std::endl;
        return 0;
    }
	// FILE *fpOut1=fopen(m_par_OceanWave.fname_WaveHeight.c_str(), "w+");
    // if(fpOut1==NULL)
    // {
    //     Info("打开文件失败: "+m_par_OceanWave.fname_WaveHeight);
    //     return 0;
    // }
	// // 向波高数据文件输入文件头
	// fprintf(fpOut1,"x[m]\ty[m]\tt[s]\th[m]\n");
    fout1<<"x[m]\ty[m]\tt[s]\th[m]\n";
	// 开始计算模拟波高数据
    printf("\nCalculating process : 1234");
	time(&start);
    // for (size_t kk = 0; kk < 3; kk++)
    // {
    //     // ui->roundProgressBar->setRange(0,100);
    //     omp_set_num_threads(8);
    //     double ind=0;
    //     #pragma omp parallel for
    //     for(int i=0;i<100;i++)
    //     {
    //         QElapsedTimer t;
    //         t.start();


    //         #pragma omp critical
    //         // cout<<i<<std::endl;;
    //         while(t.elapsed()<1000);
    //         ind++;
    //         // ui->roundProgressBar->setValue(ind);
    //     }
    // }
    
    
    // omp_set_num_threads(m_threadNumOMP);
	for(it=0;it<Nt;it++)
    {
		nProcess=(long)((it+1)*100/Nt);
		printf("\b\b\b\b%2ld%% ",nProcess);

		t=t1+it*d_t;
        ui->roundProgressBar->setRange(0,Ny);
        int ind=0;
        #pragma omp parallel for shared(omega1, d_omega, gravity,d_theta,theta1)
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
			// 			sum+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
					}	
				}
				h=sum;
			// 	// 输出计算结果
			// 	// fprintf(fpOut1,"%.3lf\t%.3lf\t%.3lf\t%.3lf\n",x,y,t,h);
            //     // fout1<<x<<" "<<y<<" "<<t<<" "<<h<<std::endl;
			}
            #pragma omp critical
            ind++;
            ui->roundProgressBar->setValue((int)ind);
		}
        
	}
	// fclose(fpOut1);
    fout1.close();
	printf("\n\nThe wave height data has been simulated!\n\n");

	// // 打开保存海水速度数据的文件
	// FILE *fpOut2=fopen(DataOut2, "w+");
	// // 向海水速度数据文件输入文件头
	// fprintf(fpOut2,"x[m]\ty[m]\tz[m]\tt[s]\tVx[m/s]\tVy[m/s]\tVz[m/s]\n");

	// // 开始计算模拟海水速度数据
	// for(it=0;it<Nt;it++)
    // {
	// 	nProcess=(long)((it+1)*100/Nt);
	// 	printf("\b\b\b\b%2ld%% ",nProcess);

	// 	t=t1+it*d_t;
	// 	for(iz=0;iz<Nz;iz++)
	// 	{
	// 		z=z1+iz*d_z;
	// 		for(iy=0;iy<Ny;iy++)
	// 		{
	// 			y=y1+iy*d_y;
	// 			for(ix=0;ix<Nx;ix++)
	// 			{
	// 				x=x1+ix*d_x;
	// 				sumx=0.0;
	// 				sumy=0.0;
	// 				sumz=0.0;
	// 				for(m=0;m<M;m++)
	// 				{
	// 					omega=omega1+m*d_omega;
	// 					k=omega*omega/gravity;
	// 					for(n=0;n<N;n++)
	// 					{
	// 						theta=theta1+n*d_theta;	
	// 						sumx+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
	// 						sumy+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*sin(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
	// 						sumz+=sqrt(2.0*PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
	// 					}	
	// 				}
	// 				Vx=sumx;
	// 				Vy=sumy;
	// 				Vz=sumz;
	// 				// 输出计算结果
	// 				fprintf(fpOut2,"%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",x,y,z,t,Vx,Vy,Vz);
	// 			}
	// 		}
	// 	}
	// }
	// fclose(fpOut2);
	// printf("\n\nThe velocity data has been simulated!\n\n");

	// delete[] num[0]; 
	// delete[] num;
	
	// time(&end);
	// printf("\n\nCalculation is finished !");
	// printf("\n\nComputing time: %.5lf seconds",difftime(end,start));
   
	// // 计算结束
    // printf("\n\nCalculating completed!\t\t\n\nPlease press any key to exit ...\n\n");
    // // getch();
    
    return 0;
}

void MainWindow::updateUILayout()
{
    QRect geometry3=ui->pushButton->geometry();
    // ui->pushButton->setGeometry(ui->groupBox_chartOptions->geometry().x()+ui->groupBox_chartOptions->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
    ui->roundProgressBar->setGeometry(geometry3);
}

void MainWindow::initParameters()
{
    // 文件路径
     ui->lineEdit->setText(QDir::currentPath()+"/波高模拟数据.txt");
     ui->lineEdit_2->setText(QDir::currentPath()+"/海水速度模拟数据.txt");
}
