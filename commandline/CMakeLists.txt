# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)


set(PROGRAM_NAME GMRSSofMPE)

project(GMRSS_MPE LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GMRSS_MPE_LIB_PATH "../Library/lib")
set(GMRSS_MPE_LIB_NAME "GMRSS_MPE")

set(GMRSS_MPE_MAJOR_VERSION 1)
set(GMRSS_MPE_MINOR_VERSION 0)
set(GMRSS_MPE_PATCH_VERSION 2)
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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#windows vs2017 doesn't need opm library link 
set(LIBomp "")
if (APPLE)
  set(LIBomp "omp")
endif()


#2. Source code dir

# add H2ONaCl source and library 
include_directories(../Library/include)
link_directories(${GMRSS_MPE_LIB_PATH} "../Library/lib" "/usr/local/lib")
aux_source_directory(. SRC_MAIN)
set( Srcs_bash ./GMRSS_MPEbash.cpp ./MultiProgressBar.cpp ./getopt.cpp)
set( Hdrs_bash ./GMRSS_MPEbash.h ./MultiProgressBar.h ./getopt.h)

# add Executable target
add_executable(${PROGRAM_NAME}  ${Srcs_bash}  main.cpp)

target_link_libraries(${PROGRAM_NAME} ${GMRSS_MPE_LIB_NAME} ${LIBomp})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../GMRSS_MPEVersion.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/GMRSS_MPEVersion.h)

message(STATUS "")
message(STATUS "GMRSS_MPE ${GMRSS_MPE_VERSION} has been configured for ${GMRSS_MPE_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
