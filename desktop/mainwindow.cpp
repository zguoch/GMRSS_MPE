#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "GMRSS_MPEVersion.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    statusBar()->showMessage(tr("G&MRSSofMPE version: ")+GMRSS_MPE_VERSION);
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

void MainWindow::on_pushButton_clicked()
{

}
