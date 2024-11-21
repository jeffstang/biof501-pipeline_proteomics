#!/bin/bash
 
# This runs our nextflow
while getopts ":rb" option;  # include both 'r' and 'b' options
do
  case $option in
  r)
    echo -e "Resuming pipeline\n"
    RESUME=TRUE
    ;;
  b)
    echo -e "Running as batch job mode, log will be output to ${LOG_FILE}\n"
    BATCH=TRUE
    ;;
  *)
    echo "Invalid option: -$OPTARG"
    exit 1
    ;;
  esac
done
 
# Then parse the run cmd
NXF_SRC_MAIN=main.nf
NEXTFLOW_CMD="nextflow run ${NXF_SRC_MAIN}"
 
if [ "$RESUME" == "TRUE" ]; then
  NEXTFLOW_CMD+=" -resume"
fi
 
if [ "$BATCH" == "TRUE" ]; then
  NEXTFLOW_CMD+=" -ansi-log false > ${LOG_FILE} 2>&1"  # Redirect output and errors to the log file
fi
 
# Then execute the run cmd
eval ${NEXTFLOW_CMD}