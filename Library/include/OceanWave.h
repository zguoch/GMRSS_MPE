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
    // 所需参数的结构体
    struct Par_OceanWave
    {
        double xmin,dx,xmax,ymin,dy,ymax,zmin,dz,zmax; //X,Y,Z方向空间范围和采样间隔(m)
        double tmin, dt, tmax; //时间方向范围和采样间隔(s)
        double minDmax_xyzt[4][3]; //xyzt四维方向的最小-增量-最大值
        double gravity;
        double alpha; //主浪方向(度)
        double omega2; //最大角频率(rad/s)
        double dOmega; //角频率采样间隔(rad/s)
        double dTheta; //方位采样间隔(度)
        double U10; //海面以上10米高度风速(m/s)
        std::string fname_WaveHeight, fname_SeawaterVelocity; //数据输出的文件名
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
            std::cout<<"===========================================\n";
        }
    };
    
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