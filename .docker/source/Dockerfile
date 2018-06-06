# rethinkroboticsopensource/sns-ik:indigo-source
# Downloads the moveit source code, install remaining debian dependencies, and builds workspace

FROM rethinkroboticsopensource/sns-ik:indigo-ci-shadow-fixed
MAINTAINER Ian McMahon git@ianthe.engineer

# Setup catkin workspace
ENV CATKIN_WS=/root/ws_sns_ik
ENV ROS_DISTRO=indigo
WORKDIR $CATKIN_WS

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

# Replacing shell with bash for later docker build commands
    mv /bin/sh /bin/sh-old && \
    ln -s /bin/bash /bin/sh

# Build repo
WORKDIR $CATKIN_WS
ENV TERM xterm
ENV PYTHONIOENCODING UTF-8
RUN catkin config --extend /opt/ros/$ROS_DISTRO --install --cmake-args -DCMAKE_BUILD_TYPE=Release && \
    catkin build
