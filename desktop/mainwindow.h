#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifdef _WIN32
  #define UNIT_T "C"
#else
  #define UNIT_T "°C"
#endif

// -------- VTK head files -------
#include <vtkContextView.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkElevationFilter.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkVectorText.h>
#include <vtkChartXY.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
// ============================
#include <iostream>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
protected:
    // vtk variable
    vtkSmartPointer<vtkContextView> m_vtkChartView;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> m_renderWindow;
public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    
    void init_Meter(Meter* meter,double min, double max, double value,int valuePrecision,
            double majorTick, double minorTick,int labelPrecision,
            QString label, QString unit,double radius=60);
    void updateMeter(Meter* meter,double value);
    void init_Meters();
    void updateMeters();
    void initRenderWindow();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
