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

#endif