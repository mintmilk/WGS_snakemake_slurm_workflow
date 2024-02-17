import os
import glob
from pathlib import Path

##############################环境变量修改#################################

# 定义基本目录和文件
INDEX_DIR = Path("index")

# 获取环境变量中的GATK路径或使用默认路径
GATK_PATH = os.environ.get('GATK_PATH', '/work/apps/tools/gatk/gatk-4.4.0.0/gatk')
GENOME = INDEX_DIR / "GRC6aGenome3032.fa"

# 使用glob找到所有以.R1.fastq.gz结尾的文件，并按名称排序
fastq_files_R1 = sorted(glob.glob("rawdata/*.R1.fastq.gz"))

# 提取样本名（假设文件名的格式为 '样本名.R1.fastq.gz'）
SAMPLES = [f.split("/")[-1].split(".R1.fastq.gz")[0] for f in fastq_files_R1]

print("所有样本名:", SAMPLES)
print("总共找到的样本名数量:", len(SAMPLES))

chromosomes = config['chromosomes']

##############################环境变量修改#################################


rule all:
    input:
        "all_final_output_SNP.filtered.vcf"

#构建索引
rule index_genome:
    input:
        fa=GENOME
    output:
        bwa_files = [str(GENOME) + x for x in [".amb", ".ann", ".bwt", ".pac", ".sa"]],
        fai = str(GENOME) + ".fai",
        dict = Path(GENOME).with_suffix(".dict")
    params:
        prefix = GENOME,
        partition = config['partitions'].get('default')
    shell:
       """
       bwa index -p {params.prefix} {input.fa}
       samtools faidx {input.fa}  
       {GATK_PATH} CreateSequenceDictionary -R {input.fa} -O {output.dict}
       """

# 使用 fastp 进行质量控制和过滤
rule fastp_qc:
    input:
        r1="rawdata/{sample}.R1.fastq.gz",
        r2="rawdata/{sample}.R2.fastq.gz"
    output:
        r1_out="fastp/{sample}_1_clean.fq.gz",
        r2_out="fastp/{sample}_2_clean.fq.gz",
        jsn="fastp/{sample}.json",
        hml="fastp/{sample}.html"
    threads: 2
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        mkdir -p fastp
        fastp -i {input.r1} -I {input.r2} -o {output.r1_out} -O {output.r2_out} -j {output.jsn} -h {output.hml} -w {threads} --length_required=50 --n_base_limit=6 --compression=6
        """

# 使用 BWA 进行序列比对
rule bwa_mem:
    input:
        ref = GENOME,
        r1 = "fastp/{sample}_1_clean.fq.gz",
        r2 = "fastp/{sample}_2_clean.fq.gz",
        amb = str(GENOME) + ".amb",
        ann = str(GENOME) + ".ann",
        bwt = str(GENOME) + ".bwt",
        pac = str(GENOME) + ".pac",
        sa = str(GENOME) + ".sa"
    output:
        bam ="bwa/{sample}.bam"
    threads: 4
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        partition = config['partitions']['bwa_mem'],
    shell:
        """
        mkdir -p bwa
        (time bwa mem -t {threads} -M -Y -R "@RG\\tID:foo_lane\\tPL:ILLUMINA\\tLB:library\\tSM:{wildcards.sample}" {input.ref} {input.r1} {input.r2} | samtools view -Sb - > {output.bam}) 2>&1 | tee {log}
        echo "** BWA MEM done **" 2>&1 | tee -a {log}
        """

# 使用 samtools 对 BAM 文件进行排序
rule sort_bam:
    input:
        bam="bwa/{sample}.bam",
        fai= str(GENOME) + ".fai"
    output:
        sorted_bam="bwa/{sample}.sorted.bam"
    threads: 4
    log:
        "logs/sort_bam/{sample}.log"
    resources:
        mem_mb=8000
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        (time samtools sort -@ {threads} -m {resources.mem_mb}M -O bam -o {output.sorted_bam} {input.bam} 2>&1; echo "** sorted raw bamfile done **") | tee {log}
        """

# 使用 GATK 进行重复序列标记
rule mark_duplicates:
    input:
        bam="bwa/{sample}.sorted.bam"
    output:
        markdup_bam="bwa/{sample}.sorted.markdup.bam",
        metrics="bwa/{sample}.markdup_metrics.txt"
    threads: 4
    log:
        "logs/mark_duplicates/{sample}.log"
    params:
        validation_stringency="LENIENT",
        max_file_handles=1000,
        partition = config['partitions'].get('default')
    shell:
        """
        (time {GATK_PATH} MarkDuplicates -I {input.bam} -M {output.metrics} -O {output.markdup_bam} --VALIDATION_STRINGENCY {params.validation_stringency} --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP {params.max_file_handles} 2>&1; echo "** {wildcards.sample}.markdup.bam done **") | tee {log}
        """

# 对标记了重复序列的 BAM 文件进行索引
rule samtools_index:
    input:
        bam="bwa/{sample}.sorted.markdup.bam"
    output:
        bai="bwa/{sample}.sorted.markdup.bam.bai"
    threads: 4
    log:
        "logs/samtools_index/{sample}.log"
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        (time samtools index {input.bam} 2>&1; echo "** {wildcards.sample}.sorted.markdup.bam index done **") | tee {log}
        """

# 使用 GATK 进行变异检测，并生成 GVCF 文件
rule haplotype_caller:
    input:
        bam="bwa/{sample}.sorted.markdup.bam",
        bai="bwa/{sample}.sorted.markdup.bam.bai",
        fai= str(GENOME) + ".fai",
        gatk_dict= Path(GENOME).with_suffix(".dict"),
        ref= GENOME
    output:
        "gatk/{sample}.HC.g.vcf.gz"
    threads:1
    log:
        "logs/haplotypecaller/{sample}.log"
    resources:
        mem_mb=16000
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        mkdir -p gatk
        (time {GATK_PATH} --java-options "-Xmx16G -Djava.io.tmpdir=./" HaplotypeCaller -R {input.ref} -I {input.bam} --emit-ref-confidence GVCF --native-pair-hmm-threads {threads} -O {output} 2>&1; echo "** {wildcards.sample}.HC.g.vcf.gz done **") | tee {log}
        """

# 生成 sample-name-map 文件
rule create_sample_name_map:
    output:
        "sample_name_map.txt"
    params:
        partition = config['partitions'].get('default')
    run:
        with open(output[0], "w") as map_file:
            for sample in SAMPLES:
                map_file.write(f"{sample}\tgatk/{sample}.HC.g.vcf.gz\n")

rule split_intervals:
    input:
        fai = str(GENOME) + ".fai"
    output:
        bed = "intervals/intervals.bed"
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        mkdir -p intervals
        awk '{{print $1 " 1 " $2}}' {input.fai} | sed -e 's/ /\\t/g' > {output.bed}
        """

rule split_bed_by_chromosome:
    input:
        bed = "intervals/intervals.bed"
    output:
        "intervals/intervals_splited/{chrom}.bed"
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        mkdir -p intervals/intervals_splited/
        grep "^{wildcards.chrom}\t" {input.bed} > {output}
        """

rule genomics_db_import:
    input:
        gvcfs = expand("gatk/{sample}.HC.g.vcf.gz", sample=SAMPLES),
        sample_name_map = "sample_name_map.txt",
        intervals = "intervals/intervals.bed"
    output:
        db=directory("genomics_db")
    threads: 2
    log:
        "logs/genomics_db_import.log"
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        {GATK_PATH} GenomicsDBImport \
            --genomicsdb-workspace-path {output.db} \
            -L {input.intervals} \
            --batch-size 50 \
            --reader-threads {threads} \
            --sample-name-map {input.sample_name_map} \
            --max-num-intervals-to-import-in-parallel 50 \
        """

# GenotypeGVCFs for each chromosome
rule genotype_gvcfs:
    input:
        db="genomics_db",
        ref=str(GENOME),
        intervals="intervals/intervals_splited/{chrom}.bed"
    output:
        vcf="vcfs/{chrom}.vcf.gz"
    log:
        "logs/genotype_gvcfs_{chrom}.log"
    params:
        partition = config['partitions'].get('default')
    shell:
        """
        {GATK_PATH} GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            -L {input.intervals}
        """

# GatherVcfs 规则
rule gather_vcf:
    input:
        vcf_files=expand("vcfs/{chrom}.vcf.gz", chrom=chromosomes)
    output:
        gzvcf="merged.vcf.gz"
    params:
        partition = config['partitions'].get('default')
    log:
        "logs/gather_vcf.log"
    shell:
        """
        {GATK_PATH} GatherVcfs \\
        $(for f in {input.vcf_files}; do echo -I $f; done) \\
        -O {output.gzvcf}
        """

# 筛选SNP、INDEL等
rule bgzip:
    input:
        gzvcf="merged.vcf.gz"
    output:
        vcf="merged.vcf"
    threads: 16
    params:
        partition = config['partitions'].get('default')
    log:
        "logs/dbgzip.log"
    shell:
        "bgzip -@ {threads} -d {input.gzvcf}"

rule select_SNP:
    input:
        vcf="merged.vcf",
        ref=GENOME,
    output:
        "all_final_output_SNP.raw.vcf"
    params:
        partition = config['partitions'].get('default')
    log:
        "logs/select_SNP.log"
    shell:
        """
        {GATK_PATH} SelectVariants -R {input.ref} -O {output} --variant {input.vcf} --select-type-to-include SNP
        echo "** selectSNP done **"
        """

rule filter_SNP:
    input:
        ref=GENOME,
        variant="all_final_output_SNP.raw.vcf"
    output:
        "all_final_output_SNP.filtered.vcf"
    params:
        partition = config['partitions'].get('default')
    log:
        "logs/filter_SNP.log"
    shell:
        """
        {GATK_PATH} VariantFiltration -R {input.ref} -O {output} --variant {input.variant} --filter-name "snp_filter" --filter-expression "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0"
        echo "** SNP filter done **"
        """