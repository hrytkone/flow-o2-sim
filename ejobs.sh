#!/bin/bash
###############################################

###############################################
# Usage and arguments
###############################################
PROG=`basename $0`
if [ $# -ne 4 ]
then
    echo "Usage: ${PROG} comment events ptDep granularity"
    exit;
fi

###############################################
# Script arguments
###############################################
export comment=$1
export noEvents=$2
export ptDep=$3
export gran=$4

###############################################
# Run settings
###############################################
export noFileToRun=10

###############################################
# Program name
###############################################
export DoWhat=toyFlow
export oname=${DoWhat}_${comment}_PtDep${ptDep}_Gran${gran}

###############################################
# Output file locations
###############################################
export Main_DIR=`pwd`
export Out_DIR=${Main_DIR}/outputs/${oname}/data
export LOG_DIR=${Main_DIR}/outputs/${oname}/logs
export OUT_ERRORS=${Main_DIR}/outputs/${oname}/errors
mkdir -p ${Out_DIR}
mkdir -p ${OUT_ERRORS}
mkdir -p ${LOG_DIR}

###############################################
# Create executible file
###############################################
cat << EOF > exec_${oname}
#!/bin/bash -f
cd ${Main_DIR}
source setup.sh
export what=${DoWhat}
export sedN=1000
export iseg=\${SLURM_ARRAY_TASK_ID}
sedN=\`expr \${sedN} + \${iseg}\`
export outfile=${Out_DIR}/${oname}_\${sedN}.root
export Log=${LOG_DIR}/${DoWhat}-\${sedN}.log
./\${what} \${outfile} ${noEvents} ${ptDep} ${gran} \${sedN} >& \${Log}
cd ${Main_DIR}
EOF
###############################################
# Make the file executable and run it in SLURM
###############################################
chmod +x exec_${oname}
## If you want to test the run, add --test-only to the beginning
sbatch -v --array=1-${noFileToRun} exec_${oname} -J ${oname}  -e ${OUT_ERRORS} -o ${OUT_ERRORS}
