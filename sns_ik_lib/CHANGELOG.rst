^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package sns_ik_lib
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0.2.3 (2017-10-29)
------------------
* CMakeLists Eigen cleanup
  In ROS Kinetic, cmake_modules is deprecated,
  so we will use some alternative CMakeLists
  strategies to find_package the Eigen 3.x library.
* Fixes Eigen scalar sum warning
  Eigen doesn't like the fact that we're creating an Array
  of Bools, and then attempting to sum those booleans up.
  Instead, we need to cast the Array into an int, and then
  sum over it to store into the sum integer.
  Resolves https://github.com/RethinkRobotics-opensource/sns_ik/issues/56
* Two small bug fixes
  1) Pass vector by reference in getJointNames
  2) Properly fill in matrix in pinv_forBarP
* Contributors: Forrest Rogers-Marcovitz, Ian McMahon

0.2.1 (2016-10-25)
------------------
* The nullspace jacobian size was transposed when using a subset of joints
* Contributors: Forrest Rogers-Marcovitz

0.2.0 (2016-09-06)
------------------
* Increase scale margin for FOSNS. Caused less unnecessary scaling with smoother motions.
* Small code cleanup
* Added NS velocity bias as task option
* Setter for dynamic loop period
* Turn on positions limits for velocity IK, but not for position IK which uses the barrier function instead.
* Setters for joint Velocity and Acceleration limits
* Correct nullspace projection for FSNS and FOSNS
* Contributors: Forrest Rogers-Marcovitz

0.1.1 (2016-04-28)
---------------------------------
* Fixed install location for sns_ik headers
* Minor code cleanup / old code removal

0.1.0 (2016-04-22)
---------------------------------
* Added setters for barrier function and added some comments
* Barrier function implemented
* Very simple limit expansion test.
* Don't divid nullspace bias by timestep in position ik solver
* Added nullspace gain to position and velocity solvers.
* Adds correct install directory and const parameters
* Turned on compile optimization -O2. Huge speed increase ~20X
* Added the standard Eigen3 package find
* Added a different way for CMake to detect eigen3 deps
* Syntax changes for stricter build warnings
* For the velocity solver in sns_ik, changed the name to CartToJntVel to minimize name confusion.
* Merge with origin/master
* First pass at nullspace bias tasks for position IK. Not optimized yet.
* First pass at nullspace bias tasks for velocity IK. Fixed P matrix for secondary tasks in fsns and fosns. Also turned off singularity logging.
* Better time comparison.
* Minor performance improvements including initializing variables outside of position for-loop and not using position limits in velocity solver when solving for position.
* Cleaner test output.
* Fixed Interface breakage & invalid function call
* Updates Position IK interface based on review comments
* Adds SharedPtr Get interface for Pos & Vel solvers
* Additional debug logging for fosns if LOG_ACTIVE is defined. Plus some additional code cleanup.
* Additional tracking of if Cartesian velocities are scaled correctly.
* Fast optimal SNS. There is still probably a bug in the code.
* Errors due to using abs instead of fabs
* Fast SNS algorithm
* Adds the ability to set velocity solver in SNS_IK
  Prior to this commit, one one need to set the velocity solver manually.
  This introduces the public setVelocitySolveType() function for creating
  and changing velocity solvers from the SNS_IK class' public interface.
  I also created some forward declarations to prevent unnecessary linking
  of velocity solver headers that aren't being used in the client's code.
* A few minor changes with parameters and logging. Random seed based on time. Small parameter changes.
* Decreased eps and switched back to standard sns vel ik
* Added testing and infra for Velocity SNS tests
* Some basic improvements to the position solver so that it stays within joint limits and stops if it gets stuck
* Added virtual symbol for proper inheritance
* Additional debug logging. Turn off later.
* Added the unused original SNS library for reference
* Adds SNS_IK() class for using URDF & limits
  This changes adds the SNS_IK class which will construct a
  KDL chain from the specified robot_description parameter.
  Position and velocity limits can be read from the URDF, and overridden
  by including a standard robot_description_planning/joint_limits yaml
  file of lower position, upper position, max velocity and max acceleration
  limits.
  Finally, a new testing executable has been included which is based largely
  off of trac_ik_examples tests. This will compare sns_ik against KDL and
  trac_ik position solvers.
* Fixed Cartesian error math which was different from Eigen
* Renamend sm class
* Disable setChain and getVelocityIK functions. Will redo with shared pointers later.
* Fixed a number of bugs in the position ik solver, though still not converging correctly
* Fixed bug in pinv_damped_P for the sub Sigma matrix
* Added Optimal SNS velocity solver with scale margin
* Correct inheritance
* Created a class for the optimal SNS velocity IK solver
* Minor change to build without warnings.
* Fixed crash related to un-initialized variables
* Attemp to add position ik test, but causing fault in velocity ik solver.
* First attempt at position IK solver
* Removed unused task elements
* Better checks for input sizes
* Merged in massive library formatting
* Converts sns lib into package framework
  This commit lays the foundation for a collection of
  SNS IK packages:
  - sns_ik : Metapackage for SNS IK solver, plugins, examples
  - sns_ik_lib: SNS IK solver library
  - sns_ik_examples: Tests and examples for SNS IK library
  - sns_ik_kinematic_plugins: A package for integrating SNS IK with moveit_ros
  Also, added Fabrezio, Ian, and Forrest as authors of these packages.
* Updated the source files to reflect sns_ik and lib move
* Moved sns_ikl to sns_ik_lib
* Contributors: Forrest Rogers-Marcovitz, Ian McMahon
