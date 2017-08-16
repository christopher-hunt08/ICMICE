#!/bin/sh
#
# multi_maus.sh
#
# Running script to automate the multi_maus submission process
#


# USER OPTIONS
################################################################################

MAUS_DIR=${MAUS_ROOT_DIR}

MULTI_MAUS="${MAUS_DIR}/files/scripts/multi_maus_skeleton.sh"

DATA_DIR="/vols/mice/ch308/dat/"

COMMAND="bin/analyze_data_offline.py"

DAQ_DATA_DIR="/vols/mice/MICE_Data/"

################################################################################


# DEFAULT OPTIONS
################################################################################

ADDITIONAL_DEFINES=""

DATA_SUBDIR=""

OUTPUT_DIR="${MAUS_DIR}/.batch_crap/"

MULTI_MAUS_CONFIG=""

DAQ_RUN_NUMBERS=""

################################################################################

if [ -z $MAUS_DIR ]
then
  echo "MAUS_ROOT_DIR not set."
  echo "Please source the environment script first."
  exit 0
fi


################################################################################
  ##  FUNCTIONS

function print_help() {
  echo
  echo "Usage : "
  echo 
  echo "multi_maus.sh <QUEUE-LENGTH> <CONFIGURATION-FILE>  [ <OPTIONS> ]"
  echo
  echo "QUEUE-LENGTH must be one of \"short\", \"medium\" or \"long\"."
  echo
  echo "Options may be any of:"
  echo 
  echo " -C     Running Command"
  echo " -O     Output Data Directory"
  echo "        (Don't use default directory name)"
  echo " -S     Output Data Subdirectory"
  echo " -D     Environment variable definition"
  echo " -R     Data run numbers"
  echo " -Q     DAQ Data Directory"
  echo
  exit 0
}


################################################################################


if (( $# )) 
then

  if [ "$1" == "-h" -o "$1" == "--help" ]
  then 
    print_help
  fi


  QUEUE_LENGTH=$1
  if [ "$QUEUE_LENGTH" != "short" ] && [ "$QUEUE_LENGTH" != "medium" ] && [ "$QUEUE_LENGTH" != "long" ]
  then
    echo "ERROR : First argument must be the chosen queue length."
    echo "Please select from \"short\", \"medium\" or \"long\""
    exit 0
  fi
  shift
else
  print_help
fi

if (( $# ))
then
  MULTI_MAUS_CONFIG=$1
  if [ ! -f $MULTI_MAUS_CONFIG ]
  then
    echo 
    echo "ERROR : Coud not find configuration file " $MULTI_MAUS_CONFIG
    echo
    exit 0
  fi
  shift
else
  echo
  echo "ERROR : Please provide the name of the configuration file"
  echo
  exit 0
fi

while [[ $# > 0 ]]
do
  key=$1

  case $key in 
    -O)
    val=$2
    if [[ $val == -* ]]
    then
      echo "ERROR : No argument supplied after option, '-O'"
      exit 0
    fi
    DATA_DIR=$val
    echo "Set data output directory to $DATA_DIR"
    shift
    shift
    ;;

    -S)
    val=$2
    if [[ $val == -* ]]
    then
      echo "ERROR : No argument supplied after option, '-S'"
      exit 0
    fi
    DATA_SUBDIR=$val
    echo "Set data subdirectory to $DATA_SUBDIR"
    shift
    shift
    ;;

    -Q)
    val=$2
    if [[ $val == -* ]]
    then
      echo "ERROR : No argument supplied after option, '-S'"
      exit 0
    fi
    DAQ_DATA_DIR=$val
    echo "Set DAQ Data directory to $DAQ_DATA_DIR"
    shift
    shift
    ;;

    -D)
    val=$2
    if [[ $val == -* ]]
    then
      echo "ERROR : No argument supplied after option, '-D'"
      exit 0
    fi
    if ! $( echo $val | egrep -q '^.*=.*$' )
    then
      echo
      echo "Invalid Variable Definition: $val"
      echo
    fi
    eval export $val
    variable_name=$(echo $val | sed 's/\(.*\)=.*/\1/g' )
    ADDITIONAL_DEFINES="$ADDITIONAL_DEFINES -v $variable_name"
    echo "Defined $val"
    shift
    shift
    ;;

    -C)
    val=$2
    if [[ $val == -* ]]
    then
      echo "ERROR : No argument supplied after option, '-C'"
      exit 0
    fi
    if [ ! -f $val ]
    then
      echo
      echo "ERROR : Could not find running command, $val"
      echo
      exit 0
    fi
    COMMAND=$val
    echo "Set command to $COMMAND"
    shift
    shift
    ;;

    -R)
    val=$2
    if [[ $val == -* ]]
    then
      echo "ERROR : No argument supplied after option, '-N'"
      exit 0
    fi
    while [ -n "$val" ] && [[ $val != -* ]]
    do
      if ! $(echo $val | egrep -q '^[0-9]{5}([[:space:]]*[0-9]{5})*$')
      then
        echo
        echo "ERROR : Run IDs must be a positive integer: '$val'"
        echo
        exit 0
      fi
      echo "Adding DAQ Run Number: $val"
      DAQ_RUN_NUMBERS="$DAQ_RUN_NUMBERS $val"
      shift
      val=$2
    done
    shift
    ;;

    *)
    echo
    echo "Warning! Unrecognised Option: $key"
    echo
    shift
    ;;
  esac

done

if [ ! -d $OUTPUT_DIR ]
then
  mkdir -p $OUTPUT_DIR
fi

DATA_FILES=""
for run_num in $DAQ_RUN_NUMBERS 
do
  FILES=`ls $DAQ_DATA_DIR/$run_num* 2> /dev/null`
  if [ -z "$FILES" ] 
  then
    echo "ERROR : $run_num Not Found"
  else
  DATA_FILES="$DATA_FILES `echo $FILES | xargs -n 1 basename | egrep '[0-9]{5}\.[0-9]{3}'`"
  fi
done

TIME=""
if [ "$QUEUE_LENGTH" == "short" ]
then 
  TIME="2:59:59"
elif [ "$QUEUE_LENGTH" == "medium" ]
then
  TIME="6:0:0"
else
  TIME="48:0:0"
fi
LENGTH="-l h_rt=$TIME"

OPTIONS="$OPTIONS $LENGTH -q hep.q -v MULTI_MAUS_MAUSDIR -v MULTI_MAUS_CONFIG -v MULTI_MAUS_COMMAND -v MULTI_MAUS_DATADIR -v MULTI_MAUS_DATA_SUBDIR -v MULTI_MAUS_DAQ_FILE -j y -o $OUTPUT_DIR $ADDITIONAL_DEFINES"

echo
echo "Using configuration file " $MULTI_MAUS_CONFIG
echo "Launching" `echo $DATA_FILES | wc -w` "Jobs to the $QUEUE queue"
echo
echo "Using the following command:"
echo
echo "qsub $OPTIONS"
echo

export MULTI_MAUS_MAUSDIR=$MAUS_DIR
export MULTI_MAUS_CONFIG
export MULTI_MAUS_COMMAND=$COMMAND
export MULTI_MAUS_DATADIR=$DATA_DIR
export MULTI_MAUS_DATA_SUBDIR=$DATA_SUBDIR

for DATA_FILE in $DATA_FILES 
do
  export MULTI_MAUS_DAQ_FILE=$DATA_FILE

#  echo "$DATA_FILE"
#  echo "qsub $OPTIONS $MULTI_MAUS"
  qsub $OPTIONS $MULTI_MAUS
done
echo
echo


################################################################################


