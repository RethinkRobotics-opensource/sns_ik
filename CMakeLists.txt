cmake_minimum_required(VERSION 2.8.3)
project(sns_ikl)

# Verify c++11 / c++0x is valid for this compiler
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "-std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++0X or C++11 support. Please choose different C++ compiler.")
endif()

find_package(catkin REQUIRED
  COMPONENTS
  roscpp
  std_msgs
  cmake_modules
)

find_package(Eigen REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  DEPENDS Eigen
  CATKIN_DEPENDS roscpp std_msgs
  LIBRARIES ${PROJECT_NAME}
)

include_directories(include ${catkin_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS})
install(DIRECTORY include/${PROJECT_NAME}/ DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})

add_library(${PROJECT_NAME} src/sns_ikl_math_utils.cpp src/sns_velocity_ik.cpp)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} ${Eigen_LIBRARIES})
install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})


add_executable(test_ikl src/test_ikl.cpp)
target_link_libraries(test_ikl ${PROJECT_NAME} ${catkin_LIBRARIES} ${Eigen_LIBRARIES})
install(TARGETS test_ikl DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})
