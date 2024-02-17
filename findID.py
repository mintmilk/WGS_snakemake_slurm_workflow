import os
import re
import glob

# 使用glob找到所有以.R1.fastq.gz结尾的文件，并按名称排序
fastq_files_R1 = sorted(glob.glob("rawdata/*.R1.fastq.gz"))

# 提取样本名（假设文件名的格式为 '样本名.R1.fastq.gz'）
SAMPLES = [f.split("/")[-1].split(".R1.fastq.gz")[0] for f in fastq_files_R1]

print("所有样本名:", SAMPLES)