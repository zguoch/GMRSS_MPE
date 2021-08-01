# Command line version

## 海浪的重磁响应模拟

* 海浪模拟
  
```bash
range_x=-50/1/50
range_y=-50/1/50
range_z=0/5/10
range_t=0/500/1000
range_xyzt=${range_x}/${range_y}/${range_z}/${range_t}
dOmega=0.01 #角频率采样间隔
omega=3 #最大角频率
U20=20 #海面以上20m的风速
alpha=45 #主浪方向
dtheta=1 #方位采样间隔
nCPU=8 #使用多少个cpu进行并行计算
outputfile=test_wave #输出文件名，只需要给定基础名称即可，后缀名会自动被添加
fmt_output=grd
GMRSSofMPE OceanWave  -R $range_xyzt -M ${dOmega}/${omega2} -U $U20 -A $alpha -T $dtheta -t $nCPU -O $outputfile -f $fmt_output
```

