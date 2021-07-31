#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "GMRSS_MPEVersion.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    statusBar()->showMessage(tr("G&MRSSofMPE version: ")+GMRSS_MPE_VERSION);
    //three meters
    init_Meters();
    updateMeters();

    //round progress bar
    // ui->roundProgressBar->setVisible(false);
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

void MainWindow::init_Meters()
{
    init_Meter(ui->meter_firstVar,5, 1000, 316,0,100,20,0, tr("Pressure"), tr("bar"));
    init_Meter(ui->meter_thirdVar,0, 100, 3.2,1,10,2,0, tr("Salinity"), tr("wt. %"));
    if(ui->comboBox->currentIndex()==0)
    {
        init_Meter(ui->meter_secondVar,0, 1000, 100,0,100,20,0, tr("Temperature"), UNIT_T);
    }else if(ui->comboBox->currentIndex()==1)
    {
        init_Meter(ui->meter_secondVar,0, 4.2, 2,2,0.5,0.1,1, tr("Enthalpy"), tr("MJ/kg"));
    }
}
void MainWindow::init_Meter(Meter* meter,double min, double max, double value,int valuePrecision,
                            double majorTick, double minorTick,int labelPrecision, QString label, QString unit,double radius)
{
    meter->setMinValue( min);
    meter->setMaxValue( max);
    meter->setValue( value );
    meter->setBackgroundColor( Qt::darkGray );
    meter->setNeedleColor( Qt::blue );
    meter->setTextColor( Qt::lightGray );
    meter->setGridColor( Qt::white );
    meter->setLabelTextColor(Qt::yellow);
    meter->setLabel( label );
    meter->setUnitsLabel( unit );
    meter->setRadius( radius );
    meter->setStartScaleAngle( 35 );
    meter->setStopScaleAngle( 325 );
    meter->setScaleStep( minorTick );
    meter->setScaleGridStep( majorTick );
    meter->setDrawValue( true );
    meter->setDrawGridValues( true );
    meter->setDrawValuePrecision( valuePrecision );
    meter->setScaleLabelPrecision( labelPrecision );
//    meter->setThresholdRange( min, threshold[0], 0 );
//    meter->setThresholdRange( threshold[0], threshold[1], 1, Qt::green );
//    meter->setThresholdRange( threshold[1], max, 2, Qt::red );
}
void MainWindow::updateMeters()
{
    QRect geo_win=this->geometry();
    QRect geo_meters=ui->group_Meters->geometry();
    QRect geo_btn=ui->pushButton_2->geometry();
    if((geo_win.width()-geo_btn.x()-geo_btn.width())<geo_meters.width())
    {
        ui->group_Meters->setVisible(false);
        return;
    }else
    {
        ui->group_Meters->setVisible(true);
        ui->group_Meters->setGeometry(geo_btn.x()+120,geo_meters.y(),geo_meters.width(),geo_meters.height());
    }

   updateMeter(ui->meter_firstVar,ui->doubleSpinBox->value());
   updateMeter(ui->meter_thirdVar,ui->doubleSpinBox_3->value()*100);
   if(ui->comboBox->currentIndex()==0)
   {
       updateMeter(ui->meter_secondVar,ui->doubleSpinBox_2->value());
   }else if(ui->comboBox->currentIndex()==1)
   {
       updateMeter(ui->meter_secondVar,ui->doubleSpinBox_2->value()/1000.0);
   }
}
void MainWindow::updateMeter(Meter* meter,double value)
{
    meter->setValue(value);
}

void MainWindow::initRenderWindow()
{
    // Geometry
    vtkNew<vtkVectorText> text;
    text->SetText("Gravity and Magnetic Response Simulation System of\nMarine Physical Environment\n(G&MRSSofMPE)\nVersion 1.0\n\nChina University of Geosciences (Wuhan)");
    vtkNew<vtkElevationFilter> elevation;
    elevation->SetInputConnection(text->GetOutputPort());
    elevation->SetLowPoint(0,0,0);
    elevation->SetHighPoint(20,0,0);

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