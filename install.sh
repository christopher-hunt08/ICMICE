#!/bin/bash
#
# Install script for the ICMICE Software package
#
# Script requires one argument, the location of an existing MAUS Install. This is REQUIRED!
#

# Setup the required variables
LOCAL_DIR=$(dirname $0)
MAUS_INSTALL=$1

echo

# Check the MAUS Install Path
if [ ! -d $MAUS_INSTALL ] || [ ! -e "$MAUS_INSTALL/README" ] || [ ! -e "$MAUS_INSTALL/env.sh" ]
then
  echo "ERROR:"
  echo "	\"$MAUS_INSTALL\" is not a valid MAUS Directory."
  echo
else

  LOG_DIR=$LOCAL_DIR/log
  CONFIG_DIR=$LOCAL_DIR/.config
  MAUS_NAME=$(basename $MAUS_INSTALL)

  echo "Found MAUS Install: $MAUS_INSTALL"
  echo "Installing within directory: $LOCAL_DIR"
  echo

  mkdir -p $LOG_DIR
  mkdir -p $CONFIG_DIR

  chmod 764 $CONFIG_DIR

  echo $MAUS_INSTALL > $CONFIG_DIR/maus

  echo "Complete"
  echo
fi

