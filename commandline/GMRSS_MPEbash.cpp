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
    if(!arg.Parse(argc, argv)) return false;
    // if(!arg.Validate()) return false;
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
    
    return true;
  }
  

}
