# rethinkroboticsopensource/sns-ik:indigo-ci
# Sets up a base image to use for running Continuous Integration on Travis

FROM ros:indigo-ros-base
MAINTAINER Ian McMahon git@ianthe.engineer

ENV TERM xterm

# Setup catkin workspace
ENV CATKIN_WS=/root/ws_sns_ik
ENV ROS_DISTRO=indigo
WORKDIR $CATKIN_WS
# Continous Integration Setting
ENV IN_DOCKER 1

# Commands are combined in single RUN statement with "apt/lists" folder removal to reduce image size
# https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#minimize-the-number-of-layers
RUN \
    mkdir src && \
    cd src && \

    # Download source so that we can get necessary dependencies
    git clone https://github.com/RethinkRobotics-opensource/sns_ik.git && \
    cd sns_ik && \
    git checkout ${ROS_DISTRO}-devel && \
    cd .. && \

    # Update apt package list as previous containers clear the cache
    apt-get -qq update && \
    apt-get -qq dist-upgrade && \

    # Install some base dependencies
    apt-get -qq install -y \
        # Some source builds require a package.xml be downloaded via wget from an external location
        wget \
        # Required for rosdep command
        sudo \
        # Preferred build tools
        python-catkin-tools \
        build-essential \
        ccache && \

    # Download all dependencies
    rosdep update && \
    rosdep install -y --from-paths . --ignore-src --rosdistro ${ROS_DISTRO} --as-root=apt:false && \

    # Remove the source code from this container. TODO: in the future we may want to keep this here for further optimization of later containers
    cd .. && \
    rm -rf src/ && \

    # Clear apt-cache to reduce image size
    rm -rf /var/lib/apt/lists/*


# Continous Integration Setting
ENV IN_DOCKER 1
