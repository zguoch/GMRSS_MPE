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
#include <string> 
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMainWindow>

#include <QThread>
#include <QtConcurrent>
#include <QFuture>
#include <QFutureWatcher>
#include <QElapsedTimer>
#include <QTranslator>

#include <thread>
#include <omp.h>

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

    QFutureWatcher<int>* watcher_;
    void busy_job_finished();
    int testjob();
    int do_busy_job();
    int m_threadNumOMP;
public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
// public slots:

//   virtual void slotOpenFile();
//   virtual void slotExit();
    void busy_job();
public:
    void initRenderWindow();

private slots:
    void on_pushButton_clicked();
    int doOceanWave();
    void updateUILayout();
private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
