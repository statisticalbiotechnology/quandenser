#!/bin/bash
export LANG=en_US.UTF-8
export LC_ALL="en_US.UTF-8"
release_dir=${HOME}/release

echo "Building Ubuntu binaries"
./manager.sh -p ubuntu -r ${release_dir}/ubuntu64 > ubuntu_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of Ubuntu binaries succeeded"
fi

echo "Building Fedora binaries"
./manager.sh -p fedora -r ${release_dir}/fedora64  > fedora_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of Fedora binaries succeeded"
fi

echo "Building CentOS binaries"
./manager.sh -p centos -r ${release_dir}/centos64  > centos_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of CentOS binaries succeeded"
fi

echo "Building native 64-bit Windows binaries"
./manager.sh -p nw64 -r ${release_dir}/win64  > nw64_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of native 64-bit Windows binaries succeeded"
fi

echo "Building native 32-bit Windows binaries"
./manager.sh -p nw32 -r ${release_dir}/win32  > nw32_output.txt 2>&1
if [[ $? -eq 0 ]]; then
  echo "Building of native 32-bit Windows binaries succeeded"
fi
