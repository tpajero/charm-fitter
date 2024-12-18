#!/bin/bash

# Run `setup-env-cvmfs.sh -h` for usage help

help() {
  echo "Setup an (LCG, compiler) view"
  echo "In case no LCG version or binary tag are specified, some defaults (with GCC) will be chosen for you."
  echo
  echo "Syntax: setup-env.sh [-h|l|b|c]"
  echo "options:"
  echo "h  Print this help."
  echo "l  LCG version to use (optional, e.g. 'LCG_105')."
  echo "b  Binary tag to use (optional, e.g. 'x86_64-el9-clang16-opt')."
  echo "c  Use Clang rather than GCC as compiler (used only if 'b' is not set)."
  echo "Example: 'source setup-env.sh -l LCG_105 -b x86_64-el9-clang16-opt'"
  echo
}

# Read the options
clang=false
lcg_setup=false
binary_setup=false
use_custom_lcg_view=false
optstring=":hl:b:c"
export OPTIND=1
while getopts ${optstring} opt; do
  case ${opt} in
    h) # Display help
      help
      return ;;
    c) # Setup Clang as compiler
      clang=true ;;
    l)
      lcg_setup=true
      export LCG_VERSION=${OPTARG} ;;
    b)
      binary_setup=true
      export BINARY_TAG=${OPTARG} ;;
    ?)
      echo "Unsupported option: ${OPTARG}"
      return ;;
  esac
done
if [[ $binary_setup = "true" && $clang = "true" ]] ; then
  echo "WARNING: The option -c is ignored as the binary tag ${BINARY_TAG} was set explicitly"
fi

# Setup the LCG version and binary tag
echo "Will build GammaCombo with the following configuration:"
if grep -q "Red Hat Enterprise Linux 9\|AlmaLinux 9" /etc/os-release; then
  echo "  OS:           RHEL9/Alma9"
  if [ ${lcg_setup} = false ]; then
    export LCG_VERSION=LCG_105
  fi
  if [ ${binary_setup} = false ]; then
    if ${clang} ; then
      export BINARY_TAG=x86_64-el9-clang16-opt
    else
      export BINARY_TAG=x86_64-el9-gcc13-opt
    fi
  fi
else
  echo "The operating system of this computer is not supported by this script."
  return 1
fi
echo "  LCG version:  ${LCG_VERSION}"
echo "  Binary tag:   ${BINARY_TAG}"
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh $LCG_VERSION $BINARY_TAG
