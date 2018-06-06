#!/bin/bash

# Software License Agreement - BSD License
#
# Inspired by MoveIt! travis https://github.com/ros-planning/moveit_core/blob/09bbc196dd4388ac8d81171620c239673b624cc4/.travis.yml
# Inspired by JSK travis https://github.com/jsk-ros-pkg/jsk_travis
# Inspired by ROS Industrial https://github.com/ros-industrial/industrial_ci
#
# Maintainer: Ian McMahon
# Author:  Dave Coleman, Isaac I. Y. Saito, Robert Haschke

# Note: ROS_REPOSITORY_PATH is no longer a valid option, use ROS_REPO. See README.md

export CI_SOURCE_PATH=$(pwd) # The repository code in this pull request that we are testing
export CI_PARENT_DIR=.travis_ci  # This is the folder name that is used in downstream repositories in order to point to this repo.
export HIT_ENDOFSCRIPT=false
export REPOSITORY_NAME=${PWD##*/}
export CATKIN_WS=/root/ws_sns_ik
echo "---"
echo "Testing branch '$TRAVIS_BRANCH' of '$REPOSITORY_NAME' on ROS '$ROS_DISTRO'"

# Helper functions
source ${CI_SOURCE_PATH}/$CI_PARENT_DIR/util.sh

# Run all CI in a Docker container
if ! [ "$IN_DOCKER" ]; then

    # Choose the correct CI container to use
    case "$ROS_REPO" in
        ros-shadow-fixed)
            export DOCKER_IMAGE=rethinkroboticsopensource/sns-ik:$ROS_DISTRO-ci-shadow-fixed
            ;;
        *)
            export DOCKER_IMAGE=rethinkroboticsopensource/sns-ik:$ROS_DISTRO-ci
            ;;
    esac
    echo "Starting Docker image: $DOCKER_IMAGE"

    # Pull first to allow us to hide console output
    docker pull $DOCKER_IMAGE > /dev/null

    # Start Docker container
    docker run \
        -e ROS_REPO \
        -e ROS_DISTRO \
        -e BEFORE_SCRIPT \
        -e CI_PARENT_DIR \
        -e CI_SOURCE_PATH \
        -e UPSTREAM_WORKSPACE \
        -e TRAVIS_BRANCH \
        -e TEST \
        -e TEST_BLACKLIST \
        -v $(pwd):/root/$REPOSITORY_NAME \
        -v $HOME/.ccache:/root/.ccache \
        $DOCKER_IMAGE \
        /bin/bash -c "cd /root/$REPOSITORY_NAME; source .travis_ci/travis.sh;"
    return_value=$?

    if [ $return_value -eq 0 ]; then
        echo "$DOCKER_IMAGE container finished successfully"
        HIT_ENDOFSCRIPT=true;
        exit 0
    fi
    echo "$DOCKER_IMAGE container finished with errors"
    exit 1 # error
fi

# If we are here, we can assume we are inside a Docker container
echo "Inside Docker container"

# Update the sources
travis_run apt-get -qq update

# Make sure the packages are up-to-date
travis_run apt-get -qq dist-upgrade

# Enable ccache
travis_run apt-get -qq install ccache
export PATH=/usr/lib/ccache:$PATH

# Install and run xvfb to allow for X11-based unittests on DISPLAY :99
travis_run apt-get -qq install xvfb mesa-utils
Xvfb -screen 0 640x480x24 :99 &
export DISPLAY=:99.0
travis_run_true glxinfo

# Setup rosdep - note: "rosdep init" is already setup in base ROS Docker image
travis_run rosdep update

# Create workspace
travis_run mkdir -p $CATKIN_WS/src
travis_run cd $CATKIN_WS/src

# Install dependencies necessary to run build using .rosinstall files
if [ ! "$UPSTREAM_WORKSPACE" ]; then
    export UPSTREAM_WORKSPACE="debian";
fi

# link in the repo we are testing
travis_run ln -s $CI_SOURCE_PATH .

# Debug: see the files in current folder
travis_run ls -a

# Run before script
if [ "${BEFORE_SCRIPT// }" != "" ]; then
    travis_run sh -c "${BEFORE_SCRIPT}";
fi

# Install source-based package dependencies
travis_run rosdep install -y -q -n --from-paths . --ignore-src --rosdistro $ROS_DISTRO

# Change to base of workspace
travis_run cd $CATKIN_WS

# Configure catkin
travis_run catkin config --extend /opt/ros/$ROS_DISTRO --install --cmake-args -DCMAKE_BUILD_TYPE=Release

# Console output fix for: "WARNING: Could not encode unicode characters"
export PYTHONIOENCODING=UTF-8

# For a command that doesn’t produce output for more than 10 minutes, prefix it with travis_run_wait
travis_run_wait 60 catkin build --no-status --summarize || exit 1

travis_run ccache -s

# Source the new built workspace
travis_run source install/setup.bash;

# Choose which packages to run tests on
echo "Test blacklist: $TEST_BLACKLIST"
echo "--------------"
TEST_PKGS=$(catkin_topological_order "$CI_SOURCE_PATH" --only-names | grep -Fvxf <(echo "$TEST_BLACKLIST" | tr ' ;,' '\n') | tr '\n' ' ')

if [ -n "$TEST_PKGS" ]; then
    TEST_PKGS="--no-deps $TEST_PKGS";
    # Fix formatting of list of packages to work correctly with Travis
    IFS=' ' read -r -a TEST_PKGS <<< "$TEST_PKGS"
fi
echo -e "Test packages: ${TEST_PKGS}"

# Run catkin package tests
travis_run catkin build --no-status --summarize --make-args tests -- ${TEST_PKGS[@]}

# Run non-catkin package tests
travis_run catkin build --catkin-make-args run_tests -- --no-status --summarize ${TEST_PKGS[@]}

# Show test results and throw error if necessary
travis_run catkin_test_results

echo "Travis script has finished successfully"
HIT_ENDOFSCRIPT=true
exit 0
