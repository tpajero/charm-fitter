#!/bin/bash

compiler="clang"

if [ "${compiler}" = "gcc" ]; then
  if grep -q "Red Hat Enterprise Linux 9" /etc/os-release; then
    export LCG_VERSION=LCG_104
    export BINARY_TAG=x86_64-el9-gcc13-opt
  elif grep -q "CentOS Linux 7" /etc/os-release; then
    export LCG_VERSION=LCG_102
    export BINARY_TAG=x86_64-centos7-gcc8-opt
  else
    echo "The operating system of this computer is not supported by this script."
    exit
  fi
elif [ "${compiler}" = "clang" ]; then
  if grep -q "Red Hat Enterprise Linux 9" /etc/os-release; then
    export LCG_VERSION=LCG_104
    export BINARY_TAG=x86_64-el9-clang16-opt
  elif grep -q "CentOS Linux 7" /etc/os-release; then
    export LCG_VERSION=LCG_102
    export BINARY_TAG=x86_64-centos7-clang12-opt
  else
    echo "The operating system of this computer is not supported by this script."
    exit
  fi
else
  echo "Compiler not found"
  return
fi

echo "  LCG version: ${LCG_VERSION}"
echo "  Binary tag: ${BINARY_TAG}"
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh $LCG_VERSION $BINARY_TAG
