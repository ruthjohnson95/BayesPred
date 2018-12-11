#!/usr/bin/env sh
#$ -cwd
#$ -j y
#$ -l h_data=5G,h_rt=5:00:00,highp
#$ -o unity_mix_sims.log
#$ -t 1-100:1

SGE_TASK_ID=1

#source /u/local/Modules/default/init/modules.sh
#module load python/2.7

for i in {1..100}
do
    COUNTER=$((COUNTER+1))
    if [[ $COUNTER -eq $SGE_TASK_ID ]]
    then
	bash simulate_hoffman.sh "1,2,3,4" $SGE_TASK_ID 
    fi
done