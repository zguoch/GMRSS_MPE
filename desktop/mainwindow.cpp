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
                    if(getPar_OceanWave(m_par_OceanWave))
                    {
                        auto future = QtConcurrent::run(this, &MainWindow::doOceanWave);
                        watcher_->setFuture(future);
                    }else
                    {
                        busy_job_finished();
                    }
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
    ui->roundProgressBar->setValue(0);
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
bool MainWindow::getPar_OceanWave(OCEANWAVE::Par_OceanWave& par)
{
    par.minDmax_xyzt[0][0]=ui->doubleSpinBox_xmin->value();
    par.minDmax_xyzt[0][1]=ui->doubleSpinBox_dx->value();
    par.minDmax_xyzt[0][2]=ui->doubleSpinBox_xmax->value();
    par.minDmax_xyzt[1][0]=ui->doubleSpinBox_ymin->value();
    par.minDmax_xyzt[1][1]=ui->doubleSpinBox_dy->value();
    par.minDmax_xyzt[1][2]=ui->doubleSpinBox_ymax->value();
    par.minDmax_xyzt[2][0]=ui->doubleSpinBox_zmin->value();
    par.minDmax_xyzt[2][1]=ui->doubleSpinBox_dz->value();
    par.minDmax_xyzt[2][2]=ui->doubleSpinBox_zmax->value();
    par.minDmax_xyzt[3][0]=ui->doubleSpinBox_tmin->value();
    par.minDmax_xyzt[3][1]=ui->doubleSpinBox_dt->value();
    par.minDmax_xyzt[3][2]=ui->doubleSpinBox_tmax->value();
    par.alpha=ui->doubleSpinBox_alpha->value();
    par.dTheta=ui->doubleSpinBox_dTheta->value();
    par.U10=ui->doubleSpinBox_U10->value();
    par.omega2=ui->doubleSpinBox_omega2->value();
    par.dOmega=ui->doubleSpinBox_dOmega->value();
    par.gravity=ui->doubleSpinBox_gravity->value();

    std::string outputPath = ui->lineEdit->text().toStdString();
    par.fname_WaveHeight=outputPath + +"/" + ui->lineEdit_2->text().toStdString() + "_h";
    par.fname_SeawaterVelocity=outputPath + +"/" + ui->lineEdit_2->text().toStdString() + "_U";
    par.fmt_outputFile = ui->comboBox_fileFormat->currentText().toStdString();
    // 测试文件路径是否正常
    std::string fname_test = outputPath+"/tmp.test";
    std::ofstream fout_test(fname_test);
    if(!fout_test)
    {
        Info("输出文件路径不正确，无法在此路径下写入文件："+outputPath);
        return false;
    }else
    {
        fout_test.close();
        std::remove(fname_test.c_str());
    }
    // -------------
    par.nThreads = m_threadNumOMP;
    return true;
}
int MainWindow::doOceanWave()
{
    OCEANWAVE::Par_OceanWave parm = m_par_OceanWave;
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
    long Nt=(long)((t2-t1)/d_t+1.0);
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
    srand(time(NULL));
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            num[i][j]=OCEANWAVE::uniform(0.0,2.0*PI);
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
    ui->roundProgressBar->setRange(0, Ny*Nt);
    int ind=0;
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
                        sum+=sqrt(2.0*OCEANWAVE::PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
                    }	
                }
                // // 根据情况判断是否需要把数据返回
                // if(returnResults == RETURN_RESULTS_LATEST && it==(Nt-1))
                // {
                //     h[0][iy][ix]=sum;
                // }else if(returnResults == RETURN_RESULTS_ALL)
                // {
                //     h[it][iy][ix]=sum;
                // }
                wave_h[iy][ix] = sum;
            }
            if(parm.showProgress){
                #pragma omp critical
                ind++;
                ui->roundProgressBar->setValue(ind);
            }
        }
        // save result at t
        OCEANWAVE::SaveResult(parm,wave_h, wave_U, t, SAVE_RESULT_h);
    }
    // 开始计算模拟海水速度数据
    ind = 0;
    ui->roundProgressBar->setRange(0, Nz*Nt);
    ui->roundProgressBar->setValue(1);
    for(it=0;it<Nt;it++)
    {
        t=t1+it*d_t;
        #pragma omp parallel for private(iz, iy, ix, m, n, omega, k, theta, z, y, x, sumx, sumy, sumz, Vx, Vy, Vz)
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
                            sumx+=sqrt(2.0*OCEANWAVE::PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
                            sumy+=sqrt(2.0*OCEANWAVE::PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*sin(theta)*sin(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
                            sumz+=sqrt(2.0*OCEANWAVE::PM2(omega,U19p5,theta-alfa,gravity)*d_omega*d_theta)*exp(-k*z)*omega*cos(k*x*cos(theta)+k*y*sin(theta)-omega*t+num[m][n]);
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
            }
            if(parm.showProgress){
                #pragma omp critical
                ind++;
                ui->roundProgressBar->setValue(ind);
            }
        }
        SaveResult(parm,wave_h, wave_U, t, SAVE_RESULT_U);
    }
    
    //release pointer
    delete[] num[0]; 
    delete[] num;

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
     ui->lineEdit->setText(QDir::currentPath());
     ui->lineEdit_2->setText("test");
}

void MainWindow::on_pushButton_2_clicked()
{
    std::string filter_ext, title_dlg;
    QString dir_output;
    title_dlg=tr("选择保存数据的路径").toStdString().c_str();
    dir_output = QFileDialog::getExistingDirectory(this, tr(title_dlg.c_str()), ui->lineEdit->text());
    // QFileDialog::getExistingDirectory();
    // QFileDialog dialog;
    // dialog.setFileMode(QFileDialog::DirectoryOnly);
    // dialog.setOption(QFileDialog::DontUseNativeDialog, true);
    // dialog.setOption(QFileDialog::ShowDirsOnly, false);
    // dialog.exec();
    if (!dir_output.isNull())
    {
        ui->lineEdit->setText(dir_output);
    }
}
