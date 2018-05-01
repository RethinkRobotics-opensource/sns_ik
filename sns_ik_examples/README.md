# README:  sns_ik_examples

This directory contains scripts that run some tests of the SNS-IK library using
the Sawyer and Baxter robot models. There are three scripts:
- `test_baxter.launch` compares the various IK solvers using the Baxter robot model
- `test_sawyer.launch` compares the various IK solvers using the Sawyer robot model
- `test_ik_solvers.launch` compares the various IK solvers on large random data sets

## Tutorial: how to run these scripts

First you need to setup your catkin workspace. If you are unfamiliar with this
process then you should refer to the ROS tutorial page:
http://wiki.ros.org/catkin/Tutorials/create_a_workspace

Once you set up the catkin workspace it should look something like:
````
ros_ws/
- src/
  - sns_ik/
  - sawyer_robot/
````

The `sns_ik/` directory contains this repository, cloned from:
https://github.com/RethinkRobotics-opensource/sns_ik.
The `sawyer_robot/` directory contains the open-source model of the Sawyer robot, cloned from:
https://github.com/RethinkRobotics/sawyer_robot.

Once the directory structure is in place you will need to setup the catkin workspace.
The commands will look something like this:
````
$ cd ros_ws
$ source /opt/ros/indigo/setup.sh
$ catkin_make
$ source devel/setup.bash
````

Now you can run the test scripts:
````
$ roslaunch sns_ik_examples test_sawyer.launch
$ roslaunch sns_ik_examples test_ik_solvers.launch
````
