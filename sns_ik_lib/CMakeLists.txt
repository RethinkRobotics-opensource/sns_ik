cmake_minimum_required(VERSION 2.8.3)
project(sns_ik_lib)

# Verify c++11 / c++0x is valid for this compiler
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -Wall")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x -O2 -Wall")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++0X or C++11 support. Please choose different C++ compiler.")
endif()

find_package(orocos_kdl REQUIRED)
find_package(Eigen3 REQUIRED)
set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIRS})
find_package(catkin REQUIRED
  COMPONENTS
  roscpp
  std_msgs
  kdl_parser
)

catkin_package(
  INCLUDE_DIRS include utilities
  DEPENDS Eigen3 orocos_kdl
  CATKIN_DEPENDS roscpp std_msgs
  LIBRARIES sns_ik
)

include_directories(utilities include ${catkin_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${orocos_kdl_INCLUDE_DIRS})
link_directories(${orocos_kdl_LIBRARY_DIRS})

# library for public API
add_library(sns_ik
            src/fosns_velocity_ik.cpp
            src/fsns_velocity_ik.cpp
            src/osns_sm_velocity_ik.cpp
            src/osns_velocity_ik.cpp
            src/sns_acc_ik_base.cpp
            src/sns_ik.cpp
            src/sns_ik_base.cpp
            src/sns_position_ik.cpp
            src/sns_vel_ik_base.cpp
            src/sns_vel_ik_base_interface.cpp
            src/sns_velocity_ik.cpp
            utilities/sns_ik_math_utils.cpp
            utilities/sns_linear_solver.cpp)
target_link_libraries(sns_ik ${catkin_LIBRARIES} ${orocos_kdl_LIBRARIES})

# install the public API
install(TARGETS sns_ik LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
install(DIRECTORY include/ DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION})
install(DIRECTORY utilities/ DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION})

# Unit tests (google framework)
if (CATKIN_ENABLE_TESTING)

  # library for test utilities
  add_library(sns_ik_test
            test/rng_utilities.cpp
            test/sawyer_model.cpp
            test/test_utilities.cpp)
  target_link_libraries(sns_ik_test sns_ik ${catkin_LIBRARIES})

  # unit tests
  catkin_add_gtest(rng_utilities_test test/rng_utilities_test.cpp)
  target_link_libraries(rng_utilities_test sns_ik sns_ik_test ${catkin_LIBRARIES})
  catkin_add_gtest(sns_ik_math_utils_test test/sns_ik_math_utils_test.cpp)
  target_link_libraries(sns_ik_math_utils_test sns_ik sns_ik_test ${catkin_LIBRARIES})
  catkin_add_gtest(sns_ik_pos_test test/sns_ik_pos_test.cpp)
  target_link_libraries(sns_ik_pos_test sns_ik sns_ik_test ${catkin_LIBRARIES})
  catkin_add_gtest(sns_ik_vel_test test/sns_ik_vel_test.cpp)
  target_link_libraries(sns_ik_vel_test sns_ik sns_ik_test ${catkin_LIBRARIES})
  catkin_add_gtest(sns_vel_ik_base_test test/sns_vel_ik_base_test.cpp)
  target_link_libraries(sns_vel_ik_base_test sns_ik sns_ik_test ${catkin_LIBRARIES})
  catkin_add_gtest(sns_acc_ik_base_test test/sns_acc_ik_base_test.cpp)
  target_link_libraries(sns_acc_ik_base_test sns_ik sns_ik_test ${catkin_LIBRARIES})

endif()
