#include "GMRSS_MPEbash.h"
// using namespace GMRSS_MPEbash;
namespace GMRSS_MPEbash
{ 
  bool bash_run(int argc, char** argv)
  {
    #ifdef _WIN32
      // set terminal as black(bg)+white(fg) model
      system("color 07"); //see https://www.geeksforgeeks.org/how-to-print-colored-text-in-c/
      GetConsoleScreenBufferInfo(m_hConsole, &csbi); 
      m_currentConsoleAttr = csbi.wAttributes;
      int width = (int)(csbi.srWindow.Right-csbi.srWindow.Left+1);
      // int height = (int)(csbi.srWindow.Bottom-csbi.srWindow.Top+1);
      if(width>119)
      {
          StartText_artASCII();
      }else
      {
          StartText();
      }
    #else
      struct winsize w;
      ioctl(0, TIOCGWINSZ, &w);
      if(w.ws_col>119)
      {
          StartText_artASCII();
      }else
      {
          StartText();
      }
    #endif
    
    // helpINFO();
    //parse arguments and check 
    cGMRSS_MPEarg arg;
    // if(!arg.Parse(argc, argv)) return false;
    // if(!arg.Validate()) return false;
    arg.runOceanWave();
    return true;
  }


  cGMRSS_MPEarg::cGMRSS_MPEarg(/* args */)
  {

  }

  cGMRSS_MPEarg::~cGMRSS_MPEarg()
  {
  }

  bool cGMRSS_MPEarg::Parse(int argc, char** argv)
  {
    return false;
    return true;
  }
  void cGMRSS_MPEarg::runOceanWave()
  {
    m_para_OceanWave.nThreads = omp_get_max_threads();
    m_para_OceanWave.gravity=9.81;
    // X方向空间范围和采样间隔(m)
    m_para_OceanWave.xmin=-50.0;
    m_para_OceanWave.xmax=+50.0;
    m_para_OceanWave.dx=1.0;
    // Y方向空间范围和采样间隔(m)
    m_para_OceanWave.ymin=-50.0;
    m_para_OceanWave.ymax=+50.0;
    m_para_OceanWave.dy=1.0;

    // Z方向空间范围和采样间隔(m)
    m_para_OceanWave.zmin=0.0;
    m_para_OceanWave.zmax=+10.0;
    m_para_OceanWave.dz=1.0;

    // 时间方向范围和采样间隔(s)
    m_para_OceanWave.tmin=0.0;
    m_para_OceanWave.tmax=+1000.0;
    m_para_OceanWave.dt=500.0;

    // 主浪方向(度)
    m_para_OceanWave.alpha=45.0;
    // 海面以上10米高度风速(m/s)
    m_para_OceanWave.U10=20.0;

    // 最大角频率(rad/s)
    m_para_OceanWave.omega2=3.0;

    // 角频率采样间隔(rad/s)
    m_para_OceanWave.dOmega=0.01;

    // 方位采样间隔(度)
    m_para_OceanWave.dTheta=1.0;
    std::vector<std::vector<std::vector<double> > > tmp_h;
    WaveHeight(m_para_OceanWave, tmp_h);
  }

}
