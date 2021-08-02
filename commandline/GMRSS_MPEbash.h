/**
 * @file GMRSS_MPEbash.h
 * @author 杜劲松，郭志馗 (jinsongdu@cug.edu.cn, zguo@geomar.de)
 * @brief 
 * @version 0.1
 * @date 2021-08-02
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef GMRSS_MPEBASH_H
#define GMRSS_MPEBASH_H
#include "GMRSS_MPEVersion.h"
#include "getopt.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

#include "OceanWave.h"
#include "MultiProgressBar.h"
#include "omp.h"

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

#ifdef _WIN32
    #include "windows.h"
    #define BLACK			0
    #define BLUE			1
    #define GREEN			2
    #define CYAN			3
    #define RED				4
    #define MAGENTA			5
    #define BROWN			6
    #define LIGHTGRAY		7
    #define DARKGRAY		8
    #define LIGHTBLUE		9
    #define LIGHTGREEN		10
    #define LIGHTCYAN		11
    #define LIGHTRED		12
    #define LIGHTMAGENTA	13
    #define YELLOW			14
    #define WHITE			15
    static HANDLE   m_hConsole=GetStdHandle(STD_OUTPUT_HANDLE);
    static WORD     m_currentConsoleAttr;
    static CONSOLE_SCREEN_BUFFER_INFO csbi;
    #define COLOR_PURPLE ""
    #define COLOR_RED "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<""
    #define COLOR_GREEN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (GREEN & 0x0F) );cout<<""
    #define COLOR_YELLOW "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<""
    #define COLOR_BLUE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BLUE & 0x0F) );cout<<""
    #define COLOR_DEFAULT "";SetConsoleTextAttribute(m_hConsole, m_currentConsoleAttr );cout<<""
    #define ERROR_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<"Error: "<<COLOR_DEFAULT
    #define WARN_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<"Warning: "<<COLOR_DEFAULT
#else
    // define color, this seems only work on MacOS and linux, doesn't work on windows
    #define ERROR_COUT "["<<"\033[31mError: "<<"\033[0m] "
    #define WARN_COUT "["<<"\033[33mWarning: "<<"\033[0m] "
    #define COLOR_PURPLE "\033[35m"
    #define COLOR_RED "\033[31m"
    #define COLOR_GREEN "\033[32m"
    #define COLOR_YELLOW "\033[33m"
    #define COLOR_BLUE "\033[34m"
    #define COLOR_DEFAULT "\033[0m"
#endif


namespace GMRSS_MPEbash
{                                                                                            
    bool bash_run(int argc, char** argv);
    
    #define MODULE_OceanWave 0
    #define MODULE_Tail 1
    #define MODULE_OceanCurrent 2

    #define MODULE_OceanWave_Name "OceanWave"
    #define MODULE_Tail_Name "Tail"
    #define MODULE_OceanCurrent_Name "OceanCurrent"
    class cGMRSS_MPEarg
    {
    private:
        bool m_havet, m_haveO, m_haveR, m_haveA, m_haveM, m_haveT, m_haveg, m_haveU, m_havef, m_haver;
        int m_ModuleIndex, m_threadNumOMP, m_valuer;
        std::string m_valueO, m_ModuleName, m_valuef;
        std::vector<std::string> m_ModuleNames, m_valueR_str, m_valueM_str;
        double m_valueA, m_valueU, m_valueg, m_valueT;
    public:
        cGMRSS_MPEarg(/* args */);
        ~cGMRSS_MPEarg();
        OCEANWAVE::Par_OceanWave m_para_OceanWave;
        bool Parse(int argc, char** argv); //Parse arguments
        bool Validate(); // validate arguments and print corresponding error information
        bool Validate_OceanWave();
        void runOceanWave();
        bool Validate_Tail();
        bool Validate_OceanCurrent();
        void StartText();
        void helpINFO(int index_module=-1);
        void helpINFO_OceanWave();
    private: 
        bool GetOptionValue(int opt, char* optarg, double& value);
        bool isNum(string str);
        int getModuleName(std::string optarg);
        vector<string> string_split (string s, string delimiter);
    };
   
    static void StartText_artASCII()
    {
        // see https://patorjk.com/software/taag/#p=display&f=ANSI%20Shadow&t=GMRSS%20of%20MPE
        std::cout<<" ██████╗ ███╗   ███╗██████╗ ███████╗███████╗     ██████╗ ███████╗    ███╗   ███╗██████╗ ███████╗\n"
                 <<"██╔════╝ ████╗ ████║██╔══██╗██╔════╝██╔════╝    ██╔═══██╗██╔════╝    ████╗ ████║██╔══██╗██╔════╝\n"
                 <<"██║  ███╗██╔████╔██║██████╔╝███████╗███████╗    ██║   ██║█████╗      ██╔████╔██║██████╔╝█████╗  \n"
                 <<"██║   ██║██║╚██╔╝██║██╔══██╗╚════██║╚════██║    ██║   ██║██╔══╝      ██║╚██╔╝██║██╔═══╝ ██╔══╝  \n"
                 <<"╚██████╔╝██║ ╚═╝ ██║██║  ██║███████║███████║    ╚██████╔╝██║         ██║ ╚═╝ ██║██║     ███████╗\n"
                 <<" ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚══════╝╚══════╝     ╚═════╝ ╚═╝         ╚═╝     ╚═╝╚═╝     ╚══════╝\n";
                                                                                                     
    }
}
#endif
