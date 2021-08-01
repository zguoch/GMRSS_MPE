/**
 * @file OceanWave.h
 * @author 杜劲松，郭志馗 (jinsongdu@cug.edu.cn, zguo@geomar.de)
 * @brief 二维海浪波高与海水流速模拟
 * @version 0.1
 * @date 2021-07-31
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef H_OCEAN_WAVE
#define H_OCEAN_WAVE

// 头文件
#include <stdio.h>
#include <math.h>
// #include <mbstring.h>
#include <stdio.h>
#include <stdlib.h>
// #include <conio.h>
#include <time.h>
#include <memory.h>
#include <string>
#include <vector>

namespace OCEANWAVE
{
    #define RETURN_RESULTS_NO 0
    #define RETURN_RESULTS_ALL 1
    #define RETURN_RESULTS_LATEST 2

    #define SAVE_RESULT_h 0
    #define SAVE_RESULT_U 1
    // 所需参数的结构体
    struct Par_OceanWave
    {
        double minDmax_xyzt[4][3]; //xyzt四维方向的最小-增量-最大值
        double gravity = 9.8;
        double alpha; //主浪方向(度)
        double omega2; //最大角频率(rad/s)
        double dOmega; //角频率采样间隔(rad/s)
        double dTheta; //方位采样间隔(度)
        double U10; //海面以上10米高度风速(m/s)
        std::string fname_WaveHeight, fname_SeawaterVelocity; //数据输出的文件名
        std::string fmt_outputFile="txt"; //输出文件格式
        int nThreads = 1;
        bool showProgress = true;
        void print()
        {
            std::cout<<"=========== 海浪模拟所使用参数 ==============\n";
            std::vector<std::string> xyzt={"x", "y", "z", "t"};
            for (size_t i = 0; i < 4; i++)
            {
                std::cout<<xyzt[i]<<"方向范围: ";
                for (size_t j = 0; j < 3; j++)
                {
                    std::cout<<minDmax_xyzt[i][j]<<" ";
                }
                std::cout<<std::endl;
            }
            std::cout<<"主浪方向: "<<alpha<<"度\n";
            std::cout<<"角频率采样间隔: "<<dOmega<<"弧度/秒\n";
            std::cout<<"最大角频率: "<<omega2<<"弧度/秒\n";
            std::cout<<"方位采样间隔: "<<dTheta<<"度\n";
            std::cout<<"海面以上10m的风速: "<<U10<<" m/s\n";
            std::cout<<"重力加速度: "<<gravity<<" m/s2\n";
            std::cout<<"波浪文件: "<<fname_WaveHeight+"."+fmt_outputFile<<"\n";
            std::cout<<"速度文件: "<<fname_SeawaterVelocity+"."+fmt_outputFile<<"\n";
            std::cout<<"使用CPU个数: "<<nThreads<<"\n";
            std::cout<<"是否显示进度条: "<<(showProgress==true ? "是":"否")<<"\n";
            std::cout<<"===========================================\n";
        }
    };
    
    // save data
    void SaveResult(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h, 
        const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U);
    // -- x,y,data格式
    void SaveResult_txt(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h,
        const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U);
    // -- rectlinear vtk grid
    void SaveResult_VTK(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h, 
        const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U);
    // -- Surfer Grid
    void SaveResult_grd(const Par_OceanWave& parm, const std::vector<std::vector<double> >& h, 
        const std::vector<std::vector<std::vector<std::vector<double> > > >& U, double t, int h_or_U);
    // 计算t0时刻的波高
    double WaveHeight(const Par_OceanWave& parm, std::vector<std::vector<std::vector<double> > >& h);
    // 服从均匀分布的随机数生成函数
    double uniform(double a, double b);
    // 服从正态分布的随机数生成函数
    double Gaussrand(double Ave, double SD);
    // P-M海浪谱函数
    double PM1(double w,double U19p5,double gravity);
    // 二维P-M海浪谱函数
    double PM2(double w,double U19p5,double theta,double gravity);
    // 10米高度风速计算任意高度风速的函数
    double Wind(double z,double U10);
};
#endif