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
    // calculation mode: 0d, 1d, 2d, 3d
    #define CALCULATION_MODE_SINGLEPOINT 0
    #define CALCULATION_MODE_ONEDIMENSION 1
    #define CALCULATION_MODE_TWODIMENSION 2
    #define CALCULATION_MODE_THREEDIMENSION 3
    
    // variable selection
    #define VARIABLE_SELECTION_PTX 0
    #define VARIABLE_SELECTION_PHX 1
    #define VARIABLE_SELECTION_T 2
    #define VARIABLE_SELECTION_P 3
    #define VARIABLE_SELECTION_X 4
    #define VARIABLE_SELECTION_H 5
    #define VARIABLE_SELECTION_PT 6
    #define VARIABLE_SELECTION_PX 7
    #define VARIABLE_SELECTION_TX 8
    #define VARIABLE_SELECTION_PH 9
    #define VARIABLE_SELECTION_HX 10

    class cGMRSS_MPEarg
    {
    private:
        
    public:
        cGMRSS_MPEarg(/* args */);
        ~cGMRSS_MPEarg();
        bool Parse(int argc, char** argv); //Parse arguments
    
    };
   
    static void StartText()
    {
        //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
        cout<<COLOR_YELLOW;       //print text in yellow color
        cout << "***************************************************\n";
        cout << "*                 program GMRSS_MPE                   *\n";
        cout << "*                 ~~~~~~~ ~~~~~~~                 *\n";
        cout << "*  Version: "<<GMRSS_MPE_VERSION<<"                     *\n";
        cout << "*                                                 *\n";
        cout << "*  Equation of state of salt-water (H2O-NaCl)     *\n";
        cout << "*  - Independent variables: PTX, PHX              *\n";
        cout << "*  - Properties: density, enthalpy, viscosity     *\n";
        cout << "*  - saturation, salinity, phase diagram          *\n";
        cout << "*  unit:                                          *\n";
        cout << "*      temperature-deg.C,        pressure-bar     *\n";
        cout << "*      salinity-wt. % NaCl,   density-kg/m3       *\n";
        cout << "*      enthalpy-kJ/kg,        viscosity-Pa s      *\n";
        cout << "*                                                 *\n";
        cout << "* (c) Zhikui Guo, GEOMAR, "<<GMRSS_MPE_DATE<<", Kiel        *\n";
        cout << "*                                                 *\n";
        cout << "***************************************************\n";
        cout << "\n";
        cout<<COLOR_DEFAULT;
                                                                                                                                                                                                                        
    }
    static void helpINFO()
    {
        string version=GMRSS_MPE_VERSION;
        string author="Zhikui Guo";
        string locus="GEOMAR, Germany";
        string email="zguo@geomar.de";
        unsigned int wordWidth=20;
        // time_t now=time(0);
        // char* now_str=ctime(&now);
        string now_str=GMRSS_MPE_DATE;

        //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
        cout<<"========================== GMRSS_MPE ==========================="<<std::endl;;
        cout<<"GMRSS_MPE, a multi-platform program for Salt-Water (H2O-NaCl) Equation of State and thermodynamic properties calculation."<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Author "<<COLOR_GREEN<<author<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Locus "<<COLOR_GREEN<<locus<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Date "<<COLOR_GREEN<<now_str<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Version "<<COLOR_GREEN<<version<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Email "<<COLOR_GREEN<<email<<COLOR_DEFAULT<<std::endl;;
        cout<<"============================================================"<<std::endl;;
        cout<<COLOR_BLUE<<"Usage: GMRSS_MPE [options]"<<COLOR_DEFAULT<<std::endl;;
        cout<<"options:"<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -h "<<COLOR_BLUE<<"List descriptions of usage and available arguments"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -v "<<COLOR_BLUE<<"Print GMRSS_MPE version number"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -D "<<COLOR_BLUE<<"Dimension: 0, 1, 2, 3. e.g.: -D2"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -V "<<COLOR_BLUE<<"Select independent variables according to -D arguments."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Combination of: T, P, X, H. e.g.: -VXT"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -T "<<COLOR_BLUE<<"Set fixed temperature value if T is not in -V option."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -P "<<COLOR_BLUE<<"Set fixed pressure value if P is not in -V option. -P316"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -X "<<COLOR_BLUE<<"Set fixed salinity value if X is not in -V option."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -H "<<COLOR_BLUE<<"Set fixed enthalpy value if H is not in -V option."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -R "<<COLOR_BLUE<<"Set range and interval of variables in -V option, must in the save order with -V option."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"e.g.: -R0/0.001/0.9/0/1/600"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -G "<<COLOR_BLUE<<"Set input filename of PTX or PHX text file for multi-points calculation"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"only used when -V0 and no -P, -X, -T or -H arguments."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"The text file with three columns, PTX or PHX are decided by -V options."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -O "<<COLOR_BLUE<<"Set out put file name, file format is determined by file extension name."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Supported file format is vtk, csv, txt."<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -n "<<COLOR_BLUE<<"If normalize the result in vtk file. Only valid when -D 3"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -t "<<COLOR_BLUE<<"Set number of thread for parallel computing."<<COLOR_DEFAULT<<std::endl;;
        cout<<"Units:"<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Temperature "<<COLOR_BLUE<<"Degree Celsius: 273.15 deg.C = 1 K (Kelvin)"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Pressure "<<COLOR_BLUE<<"bar: 1 bar = 1e5 Pa = 0.1 MPa"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Salinity "<<COLOR_BLUE<<"Weight percent in range of [0,1]: seawater is 0.032 = 3.2 wt. % NaCl"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Enthalpy "<<COLOR_BLUE<<"Specific enthalpy: kJ/kg = 1000 J/kg"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Density "<<COLOR_BLUE<<"SI: kg/m3"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Viscosity "<<COLOR_BLUE<<"SI: Pa s"<<COLOR_DEFAULT<<std::endl;;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Saturation "<<COLOR_BLUE<<"in range of [0, 1]"<<COLOR_DEFAULT<<std::endl;;
    }
    static void StartText_artASCII()
    {
        // cout<<COLOR_GREEN<<"███████╗ █████╗ ██╗  ████████╗██╗    ██╗ █████╗ ████████╗███████╗██████╗     ███████╗ ██████╗ ███████╗\n"
        // <<"██╔════╝██╔══██╗██║  ╚══██╔══╝██║    ██║██╔══██╗╚══██╔══╝██╔════╝██╔══██╗    ██╔════╝██╔═══██╗██╔════╝\n"
        // <<"███████╗███████║██║     ██║   ██║ █╗ ██║███████║   ██║   █████╗  ██████╔╝    █████╗  ██║   ██║███████╗\n"
        // <<"╚════██║██╔══██║██║     ██║   ██║███╗██║██╔══██║   ██║   ██╔══╝  ██╔══██╗    ██╔══╝  ██║   ██║╚════██║\n"
        // <<"███████║██║  ██║███████╗██║   ╚███╔███╔╝██║  ██║   ██║   ███████╗██║  ██║    ███████╗╚██████╔╝███████║\n"
        // <<"╚══════╝╚═╝  ╚═╝╚══════╝╚═╝    ╚══╝╚══╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝    ╚══════╝ ╚═════╝ ╚══════╝\n"
        // <<COLOR_DEFAULT<<std::endl;;  
        
        cout<<COLOR_GREEN<<"                              $$$$$$$$\\  $$$$$$\\   $$$$$$\\  \n"
        <<"                              $$  _____|$$  __$$\\ $$  __$$\\ \n"
        <<" $$$$$$$\\ $$\\  $$\\  $$\\       $$ |      $$ /  $$ |$$ /  \\__|\n"
        <<"$$  _____|$$ | $$ | $$ |      $$$$$\\    $$ |  $$ |\\$$$$$$\\  \n"
        <<"\\$$$$$$\\  $$ | $$ | $$ |      $$  __|   $$ |  $$ | \\____$$\\ \n"
        <<" \\____$$\\ $$ | $$ | $$ |      $$ |      $$ |  $$ |$$\\   $$ |\n"
        <<"$$$$$$$  |\\$$$$$\\$$$$  |      $$$$$$$$\\  $$$$$$  |\\$$$$$$  |\n"
        <<"\\_______/  \\_____\\____/       \\________| \\______/  \\______/ \n"
        <<COLOR_DEFAULT<<std::endl;;                                                                                                          
    }
}
#endif
