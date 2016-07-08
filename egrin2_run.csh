#!/bin/csh -f

#$ -m abe
#$ -M dreiss.isb@gmail.com.org
#$ -P Baliga
#$ -cwd
##$ -l h_rt=10000:00:00
##$ -l cores=8
#$ -t 301-312
##### Tell it you're using 2 threads, so we don't fill the whole node with 40 threads (20 cores!)
#$ -pe serial 2
#$ -o output/$TASK_ID.log 
#$ -j y

#uname -a
#pwd

##/tools/bin/R CMD BATCH "--args ID=$SGE_TASK_ID" zzz_eco5_sg.R "output/zzz_eco_${SGE_TASK_ID}.Rout"

set id = "`printf '%03d' ${SGE_TASK_ID}`"
/tools/bin/R CMD BATCH --no-save --no-restore --iter=${SGE_TASK_ID} zzz_eco5.R output/zzz_eco_${id}.Rout

