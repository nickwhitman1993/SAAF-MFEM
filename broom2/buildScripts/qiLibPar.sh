#! /bin/sh

set -e
set -o pipefail

usage ()
{
   echo
   echo Useage:
   echo "    $0 compiler_options_file"
   echo
}

if [[ $# -ne 1 ]]; then
   echo
   echo ERROR: Need to specify file with compiler options.
   usage
else
   if [ -e $1 ]; then
      source ./$1
      source ./versions.sh
      source ./QuickInstallLibOMPMPI.sh
   else
      echo
      echo ERROR: compiler options file \'$1\' not found.
      usage
   fi
fi

