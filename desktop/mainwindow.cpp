#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "GMRSS_MPEVersion.h"
// 为了对核心代码的函数进行重复使用并且能实现qt的进度条
#define OCEANWAVE_DESKTOP_APP 1

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
    // ui->splitter_3->setSizes(QList<int>() << 1 << 5);
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
            case INDEX_OCEANWAVE_GRAV:
                {
                    ui->pushButton->setVisible(true);
                    ui->roundProgressBar->setVisible(false);
                    QMessageBox msgBox;
                    msgBox.setWindowTitle(tr("Information"));
                    msgBox.setText(ui->tabWidget->tabText(index_tab)+tr("-重力模拟: 正在开发中..."));
                    msgBox.setStandardButtons(QMessageBox::Yes);
                    msgBox.setDefaultButton(QMessageBox::Yes);
                    msgBox.exec();
                }
                break;
            case INDEX_OCEANWAVE_MAG:
                {
                    ui->pushButton->setVisible(true);
                    ui->roundProgressBar->setVisible(false);
                    QMessageBox msgBox;
                    msgBox.setWindowTitle(tr("Information"));
                    msgBox.setText(ui->tabWidget->tabText(index_tab)+tr("-磁力模拟: 正在开发中..."));
                    msgBox.setStandardButtons(QMessageBox::Yes);
                    msgBox.setDefaultButton(QMessageBox::Yes);
                    msgBox.exec();
                }
                break;
            default:
                break;
            }
        }
        break;
    default:
        ui->pushButton->setVisible(true);
        ui->roundProgressBar->setVisible(false);
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
    else if(ui->radioButton_8->isChecked()) {return INDEX_OCEANWAVE_MAG;}

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
// --- 这里使用预编译对核心代码进行重复使用 ------
// 必须在mainwindow.cpp里面定义#define OCEANWAVE_DESKTOP_APP 1
#include "OceanWave_kernel.cpp"
// ==========================================
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
