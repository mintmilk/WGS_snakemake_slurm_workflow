#!/bin/bash

source /work/home/zbgroup02/miniconda3/bin/activate
echo "conda ok"

conda activate wga
echo "activate wga"

module load tools/java/v20.0.1
echo "Java-v20.0.1 ok"

export GATK_PATH="/work/apps/tools/gatk/gatk-4.4.0.0/gatk"

mkdir -p logs
snakemake -s fastp2cvf.smk --configfile config.yaml --use-conda --cluster "sbatch --partition={params.partition} --cpus-per-task={threads} --mem={resources.mem_mb}M --output=logs/slurm_{rule}.out --error=logs/slurm_{rule}.err" -j 200  --rerun-incomplete

echo "脚本执行完成" 
# curl https://api.day.app/EFoykKEWiQMTNRq7wscPQZ/CallSNP_Completed