# 0. CMake Minimum version
cmake_minimum_required (VERSION 2.12)
if (POLICY CMP0086)
  # cmake warns if loaded from a min-3.0-required parent dir, so silence the warning:
  cmake_policy(SET CMP0086 NEW)
endif()
if (POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
if (POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

project(GMRSS_MPE CXX C)
# project(${PROJECT_NAME} LANGUAGES CXX)
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

set(GMRSS_MPE_MAJOR_VERSION 1)
set(GMRSS_MPE_MINOR_VERSION 7)
set(GMRSS_MPE_PATCH_VERSION 0)
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

#2. Source code dir
aux_source_directory(. SRC_MAIN)
# Head files dir
include_directories(src include)
# source
add_subdirectory(src)

# include path 
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
# 所有模块源文件
FILE(GLOB SC_MODULES  "${PROJECT_SOURCE_DIR}/src/OceanWave.cpp")
# 将所有模块编译为一个动态库
add_library (${PROJECT_NAME}  ${SC_MODULES} )
target_link_libraries(${PROJECT_NAME})
# 安装动态库路径
install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# 添加可执行程序测试模块
add_executable(OceanWave  ${Srcs_bash}  "${PROJECT_SOURCE_DIR}/src/OceanWave_test.cpp")
target_link_libraries(OceanWave ${PROJECT_NAME})

message(STATUS "")
message(STATUS "SWEOS ${GMRSS_MPE_VERSION} has been configured for ${GMRSS_MPE_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
# ===================================================================================================


