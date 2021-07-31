cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)
foreach(p
    CMP0071 # 3.10: Let AUTOMOC and AUTOUIC process GENERATED files
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

project(GMRSS_MPE LANGUAGES CXX)
set(PROGRAM_NAME GMRSSofMPE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro()
macro(set_config_option VARNAME STRING)
  set(${VARNAME} TRUE)
  list(APPEND CONFIG_OPTIONS ${STRING})
  message(STATUS "Found " ${STRING})
endmacro()
opt(OpenMP "Enable OpenMP to parallel computing" ON)
option(USE_VTK "whether use VTK library" ON)
set(GMRSS_MPE_LIB_PATH "../Library/lib")
set(GMRSS_MPE_LIB_NAME "GMRSS_MPE")

set(GMRSS_MPE_MAJOR_VERSION 1)
set(GMRSS_MPE_MINOR_VERSION 0)
set(GMRSS_MPE_PATCH_VERSION 1)
set(GMRSS_MPE_EXTRA_VERSION "")
set(GMRSS_MPE_EXTRA_VERSION_TEXI "${GMRSS_MPE_EXTRA_VERSION}")
if(NOT GMRSS_MPE_RELEASE)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ERROR_QUIET
                    OUTPUT_VARIABLE GIT_COMMIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  if(GIT_COMMIT_HASH)
    set(GMRSS_MPE_EXTRA_VERSION "${GMRSS_MPE_EXTRA_VERSION}-git-${GIT_COMMIT_HASH}")
  endif()
  set(GMRSS_MPE_EXTRA_VERSION_TEXI "${GMRSS_MPE_EXTRA_VERSION_TEXI} (development version)")
endif()
set(GMRSS_MPE_VERSION "${GMRSS_MPE_MAJOR_VERSION}.${GMRSS_MPE_MINOR_VERSION}")
set(GMRSS_MPE_VERSION "${GMRSS_MPE_VERSION}.${GMRSS_MPE_PATCH_VERSION}${GMRSS_MPE_EXTRA_VERSION}")
set(GMRSS_MPE_SHORT_LICENSE "GNU General Public License")

if(NOT DATE)
  string(TIMESTAMP DATE "%Y-%m-%d")
endif()
set(GMRSS_MPE_DATE "${DATE}")

if(NOT GMRSS_MPE_HOST)
  execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT HOSTNAME)
    set(HOSTNAME "unknown")
  endif()
  set(GMRSS_MPE_HOST "${HOSTNAME}")
endif()

if(NOT GMRSS_MPE_PACKAGER)
  execute_process(COMMAND whoami OUTPUT_VARIABLE PACKAGER
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT PACKAGER)
    set(PACKAGER "unknown")
  endif()
  string(REPLACE "\\" " " PACKAGER ${PACKAGER})
  set(GMRSS_MPE_PACKAGER "${PACKAGER}")
endif()

if(APPLE)
  set(GMRSS_MPE_OS "MacOSX")
elseif(CYGWIN OR MSYS)
  # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
  # as a pure Windows build
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
     CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                    OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CXX_COMPILER_MACHINE MATCHES "mingw")
      set(GMRSS_MPE_OS "Windows")
      set(WIN32 1)
      add_definitions(-DWIN32)
    endif()
  endif()
else()
  set(GMRSS_MPE_OS "${CMAKE_SYSTEM_NAME}")
endif()

# Using OpenMP
include_directories(/usr/local/Cellar/libomp/11.0.0/include /usr/local/Cellar/libomp/11.1.0/include /usr/local/Cellar/libomp/12.0.1/include)
find_package(OpenMP REQUIRED)
if(APPLE)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()
endif()
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set_config_option(HAVE_OpenMP "OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#windows vs2017 doesn't need opm library link 
set(LIBomp "")
if (APPLE)
  set(LIBomp "omp")
endif()
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# set icon
if(APPLE)
  set(myApp_ICON ${CMAKE_CURRENT_SOURCE_DIR}/Icons/logo.icns)
else()
  set(myApp_ICON ${CMAKE_CURRENT_SOURCE_DIR}/icon.rc)
endif()

# VTK
if (USE_VTK)
  find_package(VTK COMPONENTS
  vtkCommonCore
  vtkFiltersCore
  vtkInfovisCore
  vtkInteractionStyle
  vtkViewsQt
  vtkChartsCore
  vtkCommonDataModel
  vtkRenderingContext2D
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkViewsContext2D
  vtkIOLegacy
  )
  message("VTK version: " "${VTK_VERSION}")
  message("VTK QT version: " "${VTK_QT_VERSION}")
  include(${VTK_USE_FILE})
  if("${VTK_VERSION}" STREQUAL "")
  message(FATAL_ERROR "Could not find VTK")
  else()
  set(GMRSS_MPE_VTK_VERSION ${VTK_VERSION})
  endif()
  if("${VTK_QT_VERSION}" STREQUAL "")
  message(FATAL_ERROR "VTK was not built with Qt")
  else()
  set(GMRSS_MPE_QT_VERSION ${VTK_QT_VERSION})
  endif()
endif()

# add H2ONaCl source and library 
include_directories(../Library/include ../commandline)
link_directories(${GMRSS_MPE_LIB_PATH} "/usr/local/lib" "${Qt5_DIR}/../../../lib")
FILE(GLOB HDRS_GMRSS_MPE_Library "${PROJECT_SOURCE_DIR}/../Library/include/*.h" "${PROJECT_SOURCE_DIR}/../Library/include/*.H")
FILE(GLOB GMRSS_MPE_LIB_FILE "${GMRSS_MPE_LIB_PATH}/*${GMRSS_MPE_LIB_NAME}*")

# Use the include path and library for Qt that is used by VTK.
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set your files and resources here
set( Srcs main.cpp  mainwindow.cpp qroundprogressbar.cpp  Meter.cpp)
set( Hdrs mainwindow.h qroundprogressbar.h  Meter.h)
set( MOC_Hdrs ui_mainwindow.h )
set( UIs mainwindow.ui )
set( QRCs Icons/icons.qrc)
set( Srcs_bash ../commandline/GMRSS_MPEbash.cpp ../commandline/MultiProgressBar.cpp ../commandline/getopt.cpp)
set( Hdrs_bash ../commandline/GMRSS_MPEbash.h ../commandline/MultiProgressBar.h ../commandline/getopt.h)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# The rest should just work (sure...)
# We have ui files, this will bring in the macro: qt5_wrap_ui
find_package(Qt5 COMPONENTS Widgets Core Concurrent LinguistTools REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})
qt5_add_resources(QRC_Srcs ${QRCs} )

#languages
set(translation_zh "zh_CN")
set(TS_FILES languages/${translation_zh}.ts)
qt5_add_translation(QON_QM_FILES ${TS_FILES})

source_group("Resources" FILES
  ${UIs}
  ${QRCs}
  ${EXE_ICON} # Not present
)
source_group("Generated" FILES
  ${UI_Srcs}
  ${MOC_Srcs}
  ${QRC_Srcs}
)

# if want to generage mac os app, add MACOSX_BUNDLE follow ${PROGRAM_NAME} 
add_executable(
  ${PROGRAM_NAME}  
  WIN32
  ${Srcs} 
  ${Hdrs} 
  ${Srcs_bash} 
  ${Hdrs_bash} 
  ${UI_Srcs} 
  ${MOC_Hdrs} 
  ${QRC_Srcs}
  ${myApp_ICON}
  ${QON_QM_FILES}
  )
# test: disable macos bundle, cancle it after debug
# set(ENABLE_OS_SPECIFIC_INSTALL ON)
if(APPLE AND NOT ENABLE_OS_SPECIFIC_INSTALL)
  message("Configure for macOS install")
  set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE ON)
  set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)
  set_source_files_properties(${QON_QM_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources/translations")
  set(MACOSX_BUNDLE_ICON_FILE logo.icns)
  set_source_files_properties(${myApp_ICON} PROPERTIES 
    MACOSX_PACKAGE_LOCATION "Resources")
endif()

target_link_libraries(${PROGRAM_NAME} Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Concurrent ${GMRSS_MPE_LIB_NAME} ${LIBFreeSteam} ${LIBomp})
if (USE_VTK)
  target_link_libraries(${PROGRAM_NAME} ${VTK_LIBRARIES} )
endif()

########################################################
###                   INSTALL                        ###
########################################################
if(WIN32)
  if(ENABLE_OS_SPECIFIC_INSTALL)
    set(GMRSS_MPE_BIN .)
    set(GMRSS_MPE_LIB .)
    set(GMRSS_MPE_DOC .)
    set(GMRSS_MPE_MAN .)
    set(GMRSS_MPE_INC .)
  else()
    include(GNUInstallDirs)
    set(GMRSS_MPE_BIN ${CMAKE_INSTALL_BINDIR})
    set(GMRSS_MPE_LIB ${CMAKE_INSTALL_LIBDIR})
    set(GMRSS_MPE_DOC ${CMAKE_INSTALL_DOCDIR})
    set(GMRSS_MPE_MAN ${CMAKE_INSTALL_MANDIR}/man1)
    set(GMRSS_MPE_INC ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
  if(CYGWIN)
    unix2dos(GMRSS_MPE_API)
    if(ENABLE_PRIVATE_API)
      unix2dos(GMRSS_MPE_PRIVATE_API)
    endif()
    unix2dos(WELCOME_FILE)
    unix2dos(SDK_FILE)
    unix2dos(LICENSE_FILE)
    unix2dos(CREDITS_FILE)
    unix2dos(CHANGELOG_FILE)
    unix2dos(TUTORIAL_FILES)
    foreach(DIR ${DEMOS_DIRS})
      file(GLOB DEMO_FILES ${DIR}/?*.*)
      unix2dos(DEMO_FILES)
    endforeach()
  endif()
elseif(APPLE AND ENABLE_OS_SPECIFIC_INSTALL)
  # set these so that the files get installed nicely in the MacOSX
  # .app bundle
  message("macOS install path config")
  set(GMRSS_MPE_BIN ../MacOS)
  set(GMRSS_MPE_LIB ../MacOS)
  set(GMRSS_MPE_DOC ../../..)
  set(GMRSS_MPE_MAN ../../..)
  set(GMRSS_MPE_INC ../MacOS)
elseif(APPLE AND NOT ENABLE_OS_SPECIFIC_INSTALL)
  # set these so that the files get installed nicely in the MacOSX
  # .app bundle
  set(GMRSS_MPE_BIN .)
  set(GMRSS_MPE_LIB ./lib)
  set(GMRSS_MPE_DOC ./doc)
  set(GMRSS_MPE_MAN ./man)
  set(GMRSS_MPE_INC ./include)
else()
  include(GNUInstallDirs)
  set(GMRSS_MPE_BIN ${CMAKE_INSTALL_BINDIR})
  set(GMRSS_MPE_LIB ${CMAKE_INSTALL_LIBDIR})
  set(GMRSS_MPE_DOC ${CMAKE_INSTALL_DOCDIR})
  set(GMRSS_MPE_MAN ${CMAKE_INSTALL_MANDIR}/man1)
  set(GMRSS_MPE_INC ${CMAKE_INSTALL_INCLUDEDIR})
endif() 

install (TARGETS ${PROGRAM_NAME} DESTINATION ${GMRSS_MPE_BIN})
# install(FILES ${Hdrs_bash} DESTINATION ${GMRSS_MPE_INC})
install(FILES ${HDRS_GMRSS_MPE_Library} DESTINATION ${GMRSS_MPE_INC})
install(FILES ${GMRSS_MPE_LIB_FILE} DESTINATION ${GMRSS_MPE_LIB})
set( LICENSE_FILE ${PROJECT_SOURCE_DIR}/../LICENSE)
set( WELCOME_FILE ${PROJECT_SOURCE_DIR}/../readme.md)
install(FILES ${LICENSE_FILE} DESTINATION .)
if(APPLE)
  install (FILES ${PROJECT_SOURCE_DIR}/../commandline/build/${PROGRAM_NAME} DESTINATION "commandline")
else()
  # install language translations
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${translation_zh}.qm DESTINATION ./Resources/translations)
  if(WIN32)
    # install vtk libs
    FILE(GLOB VTK_LIBS "${VTK_DIR}/../../../bin/*.dll")
    install (FILES ${VTK_LIBS} DESTINATION ${GMRSS_MPE_BIN})
    install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Release/ DESTINATION ${GMRSS_MPE_BIN})
    install (DIRECTORY ${PROJECT_SOURCE_DIR}/../commandline/build/Release/ DESTINATION "commandline")
  else()
    FILE(GLOB VTK_LIBS "${VTK_DIR}/../../*vtk*.*")
    install (FILES ${VTK_LIBS} DESTINATION ${GMRSS_MPE_LIB})
  endif()
  
endif()
# install API
# 0. c++
FILE(GLOB GMRSS_MPE_API_cpp "${PROJECT_SOURCE_DIR}/../test/c++/*.*")
install(FILES ${GMRSS_MPE_API_cpp} DESTINATION "API/c++/")
# 1. python
FILE(GLOB GMRSS_MPE_API_Python "${PROJECT_SOURCE_DIR}/../Library/API/python/pyGMRSS_MPE/*.*")
install(FILES ${GMRSS_MPE_API_Python} DESTINATION "API/python/")
# 2. tcl
FILE(GLOB GMRSS_MPE_API_tcl "${PROJECT_SOURCE_DIR}/../Library/API/tcl/*")
install(FILES ${GMRSS_MPE_API_tcl} DESTINATION "API/tcl/")

########################################################
###                   CPACK                          ###
########################################################
set(CPACK_PACKAGE_VENDOR "Zhikui Guo (zguo@geomar.de)")
set(CPACK_PACKAGE_VERSION_MAJOR ${GMRSS_MPE_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${GMRSS_MPE_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${GMRSS_MPE_PATCH_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_FILE ${WELCOME_FILE})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "GMRSS_MPE, a desktop application and API for salt-water (H2ONaCl) EOS and thermodynamic properties.")
if(GMRSS_MPE_EXTRA_VERSION MATCHES "-git.*") # so that we'll overwrite the archives
  set(CPACK_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${GMRSS_MPE_OS}-Installer)
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROGRAM_NAME}-git-source)
else()
  set(CPACK_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${GMRSS_MPE_VERSION}-${GMRSS_MPE_OS})
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${GMRSS_MPE_VERSION}-source)
endif()
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROGRAM_NAME}")
set(CPACK_RESOURCE_FILE_LICENSE ${LICENSE_FILE})
set(CPACK_RESOURCE_FILE_README ${WELCOME_FILE})
set(CPACK_RESOURCE_FILE_WELCOME ${WELCOME_FILE})
set(CPACK_PACKAGE_CONTACT "zguo@geomar.de")
set(CPACK_PACKAGE_EXECUTABLE "${PROGRAM_NAME}")
if(ENABLE_PACKAGE_STRIP)
  set(CPACK_STRIP_FILES TRUE)
else()
  set(CPACK_STRIP_FILES FALSE)
endif()
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}" )
if (WIN32)
    if (USE_WIX_TOOLSET)
        set(CPACK_GENERATOR "WIX") # this need WiX Tooset installed and a path to candle.exe
    elseif(USE_7Z_CPACK)
    set(CPACK_GENERATOR "7Z")
    else ()
        set(CPACK_GENERATOR "NSIS") # this needs NSIS installed, and available
    endif ()
    set(CPACK_SOURCE_GENERATOR "ZIP")
    # set icon
    set(CPACK_NSIS_DISPLAY_NAME  "${PROGRAM_NAME}")
    set(CPACK_NSIS_CONTACT "zguo@geomar.de")
    set(CPACK_PACKAGE_EXECUTABLES "${PROGRAM_NAME}" "GMRSS_MPE")
    set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${PROGRAM_NAME}-${GMRSS_MPE_VERSION}")
    set(CPACK_NSIS_MENU_LINKS "https://www.scibyte.cn" "Homepage for ${PROGRAM_NAME}")
    set(CPACK_NSIS_URL_INFO_ABOUT "https://www.scibyte.cn")
    set(CPACK_NSIS_HELP_LINK "https://www.scibyte.cn")
    set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/icon.ico")
    # set(CPACK_NSIS_MUI_HEADERIMAGE "${PROJECT_SOURCE_DIR}/../manual/home/static/images/GMRSS_MPE_mac.jpg")
    set(CPACK_NSIS_WELCOME_TITLE "GMRSS_MPE: 海洋物理环境重磁相应模拟系统")
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
elseif ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_BUNDLE_APPLE_CERT_APP "Developer ID Application: Zhikui GUO (AFCJBNSQCL)")
else ()
    set(CPACK_GENERATOR "ZIP")
    set(CPACK_SOURCE_GENERATOR "TGZ") 
endif ()

include(CPack)

# Configure versions of macro in head file

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GMRSS_MPEConfig.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/GMRSS_MPEConfig.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../GMRSS_MPEVersion.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/GMRSS_MPEVersion.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../GMRSS_MPEVersion.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/../commandline/GMRSS_MPEVersion.h)

message(STATUS "")
message(STATUS "GMRSS_MPE ${GMRSS_MPE_VERSION} has been configured for ${GMRSS_MPE_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * Qt version: " ${GMRSS_MPE_QT_VERSION})
message(STATUS " * Qt5_DIR: " ${Qt5_DIR})
if (USE_VTK)
  message(STATUS " * VTK version: " ${GMRSS_MPE_VTK_VERSION})
  message(STATUS " * VTK_DIR: " ${VTK_DIR})
endif()
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")

# solve font blurry on UI
# set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)