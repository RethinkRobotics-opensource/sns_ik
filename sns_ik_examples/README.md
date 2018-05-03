# README:  sns_ik_examples

This directory contains scripts that run some tests of the SNS-IK library using
the Sawyer and Baxter robot models.

## Tutorial: how to run these scripts

First you need to setup your catkin workspace. If you are unfamiliar with this
process then you should refer to the ROS tutorial page:
http://wiki.ros.org/catkin/Tutorials/create_a_workspace

Once you set up the catkin workspace it should look something like:
````
ros_ws/
- src/
  - sns_ik/
  - baxter_robot/
  - sawyer_robot/
````

The `sns_ik/` directory contains this repository, cloned from:
https://github.com/RethinkRobotics-opensource/sns_ik.
The `baxter_robot/` director contains the open-source model of the Baxter robot, cloned from:
https://github.com/RethinkRobotics/baxter_common.git.
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
$ roslaunch sns_ik_examples test_baxter.launch
$ roslaunch sns_ik_examples test_sawyer.launch
````

## Test architecture and parameters:

The file `test_ik_solvers.launch` is the generic launch file and is called by
both of the test scripts above. The launch file then passes a set of arguments
into the `all_ik_tests` executable. The source code for `all_ik_tests` is in
`src/ik_tests.cpp`. You can read more about ROS launch scripts here:
http://wiki.ros.org/roslaunch/XML.

## How can I run the test on a different robot?

One way to do this would be to create a `test_yourNewRobot.launch` file that is
similar to either the Sawyer or Baxter test scripts that are in this package.
