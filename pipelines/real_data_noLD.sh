#!/usr/bin/env sh
#$ -cwd
#$ -j y
#$ -l h_data=10G,h_rt=24:00:00,highp
#$ -o real_data_noLD.log

GWAS_FILE=$1 

# path to folder 
MASTER_PATH=/u/home/r/ruthjohn/ruthjohn/bayesM_noLD/unity_mixture
SCRIPT_DIR=${MASTER_PATH}/scripts 
SRC_DIR=${MASTER_PATH}/src 
DATA_DIR=${MASTER_PATH}/data 

# simulation params 
SIM_NAME=$(basename $GWAS_FILE | cut -d'.' -f1) 
TRAIT=$(basename $GWAS_FILE | cut -d'.' -f1)
BINS=100
N=100000
SEED=2018 # can replace with SGE_TASK_ID
ITS=100

DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo $DATE" Starting simulation for unity-mixture: "${SIM_NAME}

# Hoffman paths 
source /u/local/Modules/default/init/modules.sh
module load python/2.7

# data will be output to DATA_DIR 
mkdir -p $DATA_DIR 

# make histogram 
python ${SCRIPT_DIR}/plot_histogram.py --name $SIM_NAME --gwas_file $GWAS_FILE --outdir ${DATA_DIR}

# find H2 from LDSC 
#bash ${SCRIPT_DIR}/ldsc_h2.sh $GWAS_FILE $N $DATA_DIR  
LDSC_H2=$(cat ${DATA_DIR}/${TRAIT}.log | grep "Total Observed scale h2:" |  head -n 1 | cut -d':' -f2 | cut -d'(' -f1 | awk '$1=$1' )

#echo $LDSC_H2 

#python ${SRC_DIR}/mixture_em_noLD.py --name $SIM_NAME --gwas_file $GWAS_FILE --bins $BINS --N $N --seed $SEED --outdir $DATA_DIR  --its $ITS --ldsc_h2 $LDSC_H2

# make all the plots 
RESULTS_FILE=${DATA_DIR}/${SIM_NAME}.${SEED}.results 
#python ${SCRIPT_DIR}/plot_EM_histogram.py --name $SIM_NAME --results_file $RESULTS_FILE --outdir ${DATA_DIR}

