#!/usr/bin/env sh


# select steps of the pipeline 
STEPS=$1

# default is all steps 
if [ -z "$STEPS" ]
then
	STEPS="1,2,3,4"
fi 

# path to folder 
#MASTER_PATH=/u/home/r/ruthjohn/ruthjohn/bayesM_noLD/unity_mixture
MASTER_PATH=/Users/ruthiejohnson/Development/unity_mixture
SCRIPT_DIR=${MASTER_PATH}/scripts 
SRC_DIR=${MASTER_PATH}/src 
DATA_DIR=${MASTER_PATH}/data 

# simulation params 
SIM_NAME=test_identity 
P_VEC=".05,.05,.90"
BINS=3
SIGMA_G=.50 
MU_VEC="0,0"
SIGMA_VEC=".001,.1"
LD_FILE=${DATA_DIR}/identity.1000.ld 
M=1000
N=100000
SEED=2018 # can replace with SGE_TASK_ID
ITS=5

DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo $DATE" Starting simulation for unity-mixture: "${SIM_NAME}

# global paths 

# Hoffman paths 
#source /u/local/Modules/default/init/modules.sh
#module load python/2.7

# data will be output to DATA_DIR 
mkdir -p $DATA_DIR 

# STEP 1: simulate gwas 
if [[ "$STEPS" =~ "1" ]]
then
	DATE=`date '+%Y-%m-%d %H:%M:%S'`
	echo $DATE" Simulting GWAS effect sizes"
	#python ${SCRIPT_DIR}/simulate.py --name $SIM_NAME --p_vec $P_VEC --mu_vec $MU_VEC --sigma_vec $SIGMA_VEC --ld_file $LD_FILE --M $M --N $N --seed $SEED --outdir $DATA_DIR
	python ${SCRIPT_DIR}/simulate.py --name $SIM_NAME --p_vec $P_VEC --bins $BINS --ld_file $LD_FILE --M $M --N $N --seed $SEED --outdir $DATA_DIR --sigma_g $SIGMA_G 
fi


# STEP 5: run inference (with NO LD)
GWAS_FILE=${DATA_DIR}/${SIM_NAME}.${SEED}.txt 
if [[ "$STEPS" =~ "2" ]]
then
	python ${SRC_DIR}/mixture_em_noLD.py --name $SIM_NAME --gwas_file $GWAS_FILE --bins $BINS --N $N --seed $SEED --outdir $DATA_DIR  --its $ITS --ldsc_h2 $SIGMA_G 
fi


# Step 6: plot gwas effect size histogram 
if [[ "$STEPS" =~ "3" ]]
then
	python ${SCRIPT_DIR}/plot_histogram.py --name $SIM_NAME --gwas_file $GWAS_FILE --outdir ${DATA_DIR}
fi


# Step 7: plot binned effect size histogram 
if [[ "$STEPS" =~ "4" ]]
then
	RESULTS_FILE=${DATA_DIR}/${SIM_NAME}.${SEED}.results 
	python ${SCRIPT_DIR}/plot_EM_histogram.py --name $SIM_NAME --results_file $RESULTS_FILE --outdir ${DATA_DIR}
fi

