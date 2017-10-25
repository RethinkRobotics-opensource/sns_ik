cmake_minimum_required(VERSION 2.8.3)
project(sns_ik_examples)

# Verify c++11 / c++0x is valid for this compiler
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "-std=c++11 -O2")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x -O2")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++0X or C++11 support. Please choose different C++ compiler.")
endif()

find_package(catkin REQUIRED
  COMPONENTS
  roscpp
  sns_ik_lib
  #trac_ik_lib #Uncomment this to test against trac_ik
)

find_package(Eigen3 REQUIRED)
set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIRS})

catkin_package(
  DEPENDS Eigen3
  CATKIN_DEPENDS roscpp
)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS}
)

add_executable(all_ik_tests src/ik_tests.cpp)
target_link_libraries(all_ik_tests ${catkin_LIBRARIES})
install(TARGETS all_ik_tests DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})
