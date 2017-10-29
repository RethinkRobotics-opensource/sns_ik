^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package sns_ik_examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0.2.3 (2017-10-29)
------------------
* CMakeLists Eigen cleanup
  In ROS Kinetic, cmake_modules is deprecated,
  so we will use some alternative CMakeLists
  strategies to find_package the Eigen 3.x library.
* Contributors: Ian McMahon

0.2.1 (2016-10-25)
------------------
* No Updates

0.2.0 (2016-09-06)
------------------
* Fixes for stricter build rules
* Contributors: Forrest Rogers-Marcovitz

0.1.1 (2016-04-28)
---------------------------------
* No updates

0.1.0 (2016-04-22)
---------------------------------
* Fixed nullspace evalution bugs
* Updated per review comments for Nullspace Tests
* Added nullspace influence calculation
* Added position limit checks on IK solvers
* Added final pose checks for position solvers
* Nullspace testing wip
* Basic nullspace bias task functionality
* Output time std deviation with plus/minus sign. Increased number of position tests.
* Turned on compile optimization -O2. Huge speed increase ~20X
* For the velocity solver in sns_ik, changed the name to CartToJntVel to minimize name confusion.
* Merge with origin/master
* First pass at nullspace bias tasks for position IK. Not optimized yet.
* Added Pos Solver Time Standard Deviation
  This commit adds Standard Deviation calculation for
  the time each position solver takes to run, for KDL,
  Trac & all SNS solvers.
* Added Launch file heierchy for multiple robots
  Now the test_ik_solvers.launch is expecting a URDF to
  have already been loaded into /robot_description or
  the specfied urdf_param location.
* Standardized seed-testing interface
* First pass at nullspace bias tasks for velocity IK. Fixed P matrix for secondary tasks in fsns and fosns. Also turned off singularity logging.
* Updated based on review comments
* Adds close delta seeds for Position IK Testing
  Adds a seed that is close to the random joint
  angle used in constructing an inverse kinematic solition.
  This should simulate the most favorable (and realistic)
  conditions for using the solver in cartesian movement.
  TBD: - Parametrize the 0.2 radian delta for the seed
  - Test how close the resulting joint solutions are to the delta
  seed
* Better time comparison.
* Allow random seeds and logging cleanup
* Check rotational speed scaling
* Cleaner test output.
* Fixed Interface breakage & invalid function call
* Updates Position IK interface based on review comments
* Adds SharedPtr Get interface for Pos & Vel solvers
* Additional tracking of if Cartesian velocities are scaled correctly.
* Fast optimal SNS. There is still probably a bug in the code.
* Fast SNS algorithm
* A few minor changes with parameters and logging. Random seed based on time. Small parameter changes.
* Added testing and infra for Velocity SNS tests
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
* Merged from master
* Fixed a number of bugs in the position ik solver, though still not converging correctly
* Added Optimal SNS velocity solver with scale margin
* Correct inheritance
* Created a class for the optimal SNS velocity IK solver
* Fixed crash related to un-initialized variables
* Attemp to add position ik test, but causing fault in velocity ik solver.
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
* Contributors: Forrest Rogers-Marcovitz, Ian McMahon
