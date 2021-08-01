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
#include <QMessageBox>
#include <QThread>
#include <QtConcurrent>
#include <QFuture>
#include <QFutureWatcher>
#include <QElapsedTimer>
#include <QTranslator>
#include <QFileDialog>

#include <thread>
#include <omp.h>

// 模拟计算库的头文件
#include "proj.h"
#include "OceanWave.h"

// INDEX define
#define INDEX_GROUP_OCEANWAVE 0
#define INDEX_GROUP_TRAIL 1

#define INDEX_OCEANWAVE_WAVE 0
#define INDEX_OCEANWAVE_GRAV 1
#define INDEX_OCEANWAVE_MAG 2

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
    OCEANWAVE::Par_OceanWave m_par_OceanWave;
public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
// public slots:

//   virtual void slotOpenFile();
//   virtual void slotExit();
    void busy_job();
    void Info(std::string infoText);
public:
    void initRenderWindow();

private slots:
    void on_pushButton_clicked();
    // 1.0 海浪模拟
    int getIndex_OceanWave_Grav_Mag();
    bool getPar_OceanWave(OCEANWAVE::Par_OceanWave& par);
    int doOceanWave();

    // 1.1 海浪的重力响应
    // ...
    // ...
    void updateUILayout();
    void initParameters();
    void on_pushButton_2_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
