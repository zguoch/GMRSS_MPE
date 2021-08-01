#include "GMRSS_MPEbash.h"
// using namespace GMRSS_MPEbash;
namespace GMRSS_MPEbash
{ 
  bool bash_run(int argc, char** argv)
  {
    if(argc==1)
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
      helpINFO();
    }
    //parse arguments and check 
    cGMRSS_MPEarg arg;
    if(!arg.Parse(argc, argv)) return false;
    if(!arg.Validate()) return false;
    return true;
  }


  cGMRSS_MPEarg::cGMRSS_MPEarg(/* args */)
  {
    m_ModuleNames.resize(5);//暂时先用一个较大的数字
    m_ModuleNames[MODULE_OceanWave]=MODULE_OceanWave_Name; //海浪
    m_ModuleNames[MODULE_Tail]=MODULE_Tail_Name; //尾迹
    m_ModuleNames[MODULE_OceanCurrent]=MODULE_OceanCurrent_Name; //洋流
  }

  cGMRSS_MPEarg::~cGMRSS_MPEarg()
  {
  }

  bool cGMRSS_MPEarg::Parse(int argc, char** argv)
  {
    if(argc<2)return false; //there is no arguments
    // check and get module name
    m_ModuleIndex = getModuleName(argv[1]);
    if(m_ModuleIndex<0)
    {
      helpINFO();
      cout<<ERROR_COUT<<"必须指定模块名: ["<<m_ModuleNames[0];
      for (size_t i = 0; i < m_ModuleNames.size(); i++)
      {
        std::cout<<"|"<<m_ModuleNames[i];
      }std::cout<<"]\n";
      
      exit(0);
    }else
    {
      const char *have_module = "-m";
      argv[1] = const_cast<char*>(have_module);//place holder, nothing to be used
    }
    int opt; 
    const char *optstring = "f:U:g:T:A:M:R:O:t:vhm"; // set argument templete
    int option_index = 0;
    int valid_args=0;
    double doubleOptValue;
    while ((opt = getopt(argc, argv, optstring)) != -1) 
    {
      if(opt!='?')
      {
        valid_args++;
      }else
      {
        cout<<WARN_COUT<<"Unknown option "<<argv[optind-1]<<endl;
      }
      switch (opt)
      {
      case 'h':
        helpINFO();
        exit(0);
        break;
      case 'v':
        cout<<"Version: "<<VERSION_MAJOR<<"."<<VERSION_MINOR<<endl;
        exit(0);
        break;
      case 't':
        m_havet=true;
        if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
        m_threadNumOMP=(int)doubleOptValue;
        if(m_threadNumOMP>omp_get_max_threads())m_threadNumOMP=omp_get_max_threads();
        if(m_threadNumOMP<1)m_threadNumOMP=1;
        break;
      case 'O':
        m_haveO=true;
        m_valueO=optarg;
        break;
      case 'R':
        m_haveR=true;
        m_valueR_str= string_split(optarg,"/");
        break;
      case 'M':   //dOmega/omega2 in OceanWave simulation
        m_haveM=true;
        m_valueM_str= string_split(optarg,"/");
        break;
      case 'A': //alpha（主浪方向) in OceanWave simulation
        m_haveA=true;
        if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
        m_valueA=doubleOptValue;
        break;
      case 'T': //dTheta in OceanWave simulation
        m_haveT=true;
        if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
        m_valueT=doubleOptValue;
        break;
      case 'g': //gravity
        m_haveg=true;
        if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
        m_valueg=doubleOptValue;
        break;
      case 'U': //U10
        m_haveU=true;
        if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
        m_valueU=doubleOptValue;
        break;
      case 'f': //file format
        m_havef=true;
        m_valuef=optarg;
        break;
      default:
        break;
      }
    }
    return true;
  }
  int cGMRSS_MPEarg::getModuleName(std::string optarg)
  {
    for (size_t i = 0; i < m_ModuleNames.size(); i++)
    {
      if(optarg==m_ModuleNames[i])
      {
        return i;
      }
    }
    return -1;
  }
  bool cGMRSS_MPEarg::Validate_OceanWave()
  {
    if(m_havef)m_para_OceanWave.fmt_outputFile = m_valuef;
    if (!m_haveO)
    {
      cout<<WARN_COUT<<"没有指定输出文件名，使用默认文件名(OceanWave_U和OceanWave_h)且保存在当前程序的运行路径"<<endl;
      m_para_OceanWave.fname_SeawaterVelocity = "OceanWave_U";
      m_para_OceanWave.fname_WaveHeight = "OceanWave_h";
    }else
    {
      m_para_OceanWave.fname_SeawaterVelocity = m_valueO+"_U";
      m_para_OceanWave.fname_WaveHeight = m_valueO+"_h";
    }
    // 检查x,y,z,t范围
    if(m_valueR_str.size()%3 != 0 || m_valueR_str.size()!=12)
    {
      cout<<ERROR_COUT<<"Option of -R argument must be a multiple of 3 and <=12, in format of [min/delta/max]"<<endl;
      return false;
    }else 
    {
      for (size_t i = 0; i < m_valueR_str.size(); i++)
      {
        if(isNum(m_valueR_str[i]))
        {
          m_para_OceanWave.minDmax_xyzt[i/3][i%3]=atof(m_valueR_str[i].c_str());
        }else
        {
          cout<<ERROR_COUT<<"The "<<i+1<<"th value in -R option is not a number: "<<COLOR_RED<<m_valueR_str[i]<<COLOR_DEFAULT<<endl;
          return false;
        }
      }
      // cout<<"计算范围\n";
      // for (size_t i = 0; i < 4; i++)
      // {
      //   for (size_t j = 0; j < 3; j++)
      //   {
      //     cout<<m_para_OceanWave.minDmax_xyzt[i][j]<<" ";
      //   }
      //   cout<<endl;
      // }
      
    }
    // 检查dOmega/omega2
    if(m_valueM_str.size()!=2)
    {
      cout<<ERROR_COUT<<"请通过-M给定dOmega(角频率采样间隔)/omega2(最大角频率), e.g. -M 0.01/3"<<endl;
      return false;
    }else
    {
      m_para_OceanWave.dOmega=atof(m_valueM_str[0].c_str());
      m_para_OceanWave.omega2=atof(m_valueM_str[1].c_str());
    }
    // 检查U10
    if(!m_haveU)
    {
      cout<<ERROR_COUT<<"请通过-U给定海面以上10m的风速, e.g. -U 20"<<endl;
      return false;
    }else
    {
      m_para_OceanWave.U10 = m_valueU;
    }
    //检查方位采样间隔
    if(!m_haveT)
    {
      cout<<ERROR_COUT<<"请通过-T给定方位采样间隔, e.g. -T 1"<<endl;
      return false;
    }else
    {
      m_para_OceanWave.dTheta = m_valueT;
    }
    //重力加速度
    if (m_haveg)m_para_OceanWave.gravity = m_valueg;
    if(m_havet)m_para_OceanWave.nThreads = m_threadNumOMP;
    //检查主浪方向
    if(!m_haveA)
    {
      cout<<ERROR_COUT<<"请通过-A给定主浪方向alpha, e.g. -A 45"<<endl;
      return false;
    }
    {
      m_para_OceanWave.alpha = m_valueA;
    }
    m_para_OceanWave.print();
    // 参数没问题那就运行吧
    std::vector<std::vector<std::vector<double> > > tmp_h;
    WaveHeight(m_para_OceanWave, tmp_h); 
    // runOceanWave();
    return true;
  }
  bool cGMRSS_MPEarg::Validate_Tail()
  {
    return true;
  }
  bool cGMRSS_MPEarg::Validate_OceanCurrent()
  {
    return true;
  }
  bool cGMRSS_MPEarg::Validate()
  {
    // 根据不同的模块调用不同的检查函数
    switch (m_ModuleIndex)
    {
    case MODULE_OceanWave:
      Validate_OceanWave();
      break;
    case MODULE_Tail:
      Validate_Tail();
      break;
    case MODULE_OceanCurrent:
      Validate_OceanCurrent();
      break;
    default:
      break;
    }
    
    return true;
  }

  bool cGMRSS_MPEarg::GetOptionValue(int opt, char* optarg, double& value)
  {
    string optarg_str=optarg;
    if(isNum(optarg_str))
    {
      value=atof(optarg);
    }else
    { 
      char optCh=opt;
      cout<<ERROR_COUT<<"Option of -"<<optCh<<" argument is empty or cannot be recognized"<<endl;
      return false;
    }
    return true;
  }
  bool cGMRSS_MPEarg::isNum(string str)
  {
      stringstream sin(str);
      double d;
      char c;
      if(!(sin >> d))
          return false;
      if (sin >> c)
          return false;
      return true;
  }
  void cGMRSS_MPEarg::runOceanWave()
  {
    m_para_OceanWave.nThreads = omp_get_max_threads();
    m_para_OceanWave.gravity=9.81;
    // X方向空间范围和采样间隔(m)
    m_para_OceanWave.minDmax_xyzt[0][0]=-50.0;
    m_para_OceanWave.minDmax_xyzt[0][2]=+50.0;
    m_para_OceanWave.minDmax_xyzt[0][1]=1.0;
    // Y方向空间范围和采样间隔(m)
    m_para_OceanWave.minDmax_xyzt[1][0]=-50.0;
    m_para_OceanWave.minDmax_xyzt[1][2]=+50.0;
    m_para_OceanWave.minDmax_xyzt[1][1]=1.0;

    // Z方向空间范围和采样间隔(m)
    m_para_OceanWave.minDmax_xyzt[2][0]=0.0;
    m_para_OceanWave.minDmax_xyzt[2][2]=+10.0;
    m_para_OceanWave.minDmax_xyzt[2][1]=1.0;

    // 时间方向范围和采样间隔(s)
    m_para_OceanWave.minDmax_xyzt[3][0]=0.0;
    m_para_OceanWave.minDmax_xyzt[3][2]=+1000.0;
    m_para_OceanWave.minDmax_xyzt[3][1]=500.0;

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

  vector<string> cGMRSS_MPEarg::string_split (string s, string delimiter) 
  {
      size_t pos_start = 0, pos_end, delim_len = delimiter.length();
      string token;
      vector<string> res;
      while ((pos_end = s.find (delimiter, pos_start)) != string::npos) 
      {
          token = s.substr (pos_start, pos_end - pos_start);
          pos_start = pos_end + delim_len;
          res.push_back (token);
      }
      res.push_back (s.substr (pos_start));
      return res;
  }
}
