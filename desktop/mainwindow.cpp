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
    auto future = QtConcurrent::run(this, &MainWindow::doOceanWave);
    watcher_->setFuture(future);
    // switch (m_dimension) {
    //     case 1:
    //     {
    //         ui->roundProgressBar->setVisible(true);
    //         auto future = QtConcurrent::run(this, &MainWindow::Calculate_Diagram1D);
    //         watcher_->setFuture(future);
    //     }
    //     break;
    //     case 2:
    //     {
    //         ui->roundProgressBar->setVisible(true);
    //         auto future = QtConcurrent::run(this, &MainWindow::Calculate_Diagram2D);
    //         watcher_->setFuture(future);
    //     }
    //     break;
    // case 3:
    //     QMessageBox msgBox;
    //     msgBox.setWindowTitle(tr("Information"));
    //     msgBox.setText(tr("3D is comming soon"));
    //     msgBox.setStandardButtons(QMessageBox::Yes);
    //     msgBox.setDefaultButton(QMessageBox::Yes);
    //     msgBox.exec();
    //     //ui->textEdit->append("3D is comming soon");
    //     break;
    // }
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

int MainWindow::doOceanWave()
{
    ui->roundProgressBar->setRange(0,500);
    int ind=0;
    omp_set_num_threads(m_threadNumOMP);
    #pragma omp parallel for
    for (size_t i = 0; i < 500; i++)
    {
        for (size_t j = 0; j < 999999; j++)
        {
            double a=999*999;
        }
        #pragma omp critical
        ind++;
        ui->roundProgressBar->setValue((int)ind);
    }
    
    return 0;
}

void MainWindow::updateUILayout()
{
    QRect geometry3=ui->pushButton->geometry();
    // ui->pushButton->setGeometry(ui->groupBox_chartOptions->geometry().x()+ui->groupBox_chartOptions->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
    ui->roundProgressBar->setGeometry(geometry3);
}