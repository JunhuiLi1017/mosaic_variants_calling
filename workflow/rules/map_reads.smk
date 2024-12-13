rule map_raw:
    input:
        get_clean_fastq
    output:
        o1=temp("{outpath}/{sample}/02_map/bwa/{sample}.raw.bam")
    log:
        l1="{outpath}/{sample}/02_map/logs/{sample}.bwa.log"
    threads:
        8
    params:
        mem="8G",
        rg=get_reads_group,
        ref=config['reference']
    shell:
        """
        bwa mem -t {threads} -M {params.rg} {params.ref} {input[0]} {input[1]} | samtools view -b -o {output.o1} > {log.l1} 2>&1
        """

rule map_raw_sort:
    input:
        "{outpath}/{sample}/02_map/bwa/{sample}.raw.bam"
    output:
        sort_bam="{outpath}/{sample}/02_map/bwa/{sample}.sort.bam",
        sort_stat="{outpath}/{sample}/02_map/bwa/{sample}.sort.stat"
    log:
        "{outpath}/{sample}/02_map/logs/{sample}.raw.sort.log"
    threads:
        8
    shell:
        """
        samtools sort -@ {threads} -m 8G -O bam -o {output.sort_bam} {input}> {log} 2>&1
        samtools stats {output.sort_bam} > {output.sort_stat}
        samtools index {output.sort_bam}
        """

rule remove_dup:
    input:
        "{outpath}/{sample}/02_map/bwa/{sample}.sort.bam"
    output:
        o1="{outpath}/02_map/dup/{sample}.sort.rmdup.bam",
        o2="{outpath}/02_map/dup/{sample}.sort.rmdup.matrix"
    log:
        "{outpath}/logs/bwa/{sample}.sort.rmdup.log"
    threads:
        8
    params:
        dedup="--REMOVE_DUPLICATES true" if config["pcr_based"] == True else "--REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
    shell:
        """
        picard MarkDuplicates --INPUT {input} --METRICS_FILE {output.o2} --OUTPUT {output.o1} {params.dedup} && samtools index {output.o1} > {log} 2>&1
        """

rule BaseRecalibrator:
    input:
        "{outpath}/{sample}/02_map/dup/{sample}.sort.rmdup.bam"
    output:
        o1="{outpath}/{sample}/02_map/bqsr/{sample}.recal_data.table"
    log:
        "{outpath}/logs/bwa/{sample}.BQSR.log"
    threads:
        8
    params:
        mem="4000",
        ref=config['reference'],
        dpsnp138=config['dpsnp138'],
        known_indels=config['known_indels'],
        Mills_and_1000G=config['Mills_and_1000G']
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar /home/junhui.li11-umw/anaconda3/envs/snakemake/share/gatk4-4.1.8.1-0/gatk-package-4.1.8.1-local.jar BaseRecalibrator -I {input} -O {output.o1} -R {params.ref} --known-sites {params.dpsnp138} --known-sites {params.known_indels} --known-sites {params.Mills_and_1000G} > {log} 2>&1
        """

rule ApplyBQSR:
    input:
        i1="{outpath}/{sample}/02_map/dup/{sample}.sort.rmdup.bam",
        i2="{outpath}/{sample}/02_map/bqsr/{sample}.recal_data.table"
    output:
        o1="{outpath}/{sample}/02_map/bqsr/{sample}.sort.rmdup.bqsr.bam",
        o2="{outpath}/{sample}/02_map/bqsr/{sample}.sort.rmdup.bqsr.bai"
    log:
        "{outpath}/logs/bwa/{sample}.applyBQSR.log"
    threads:
        8
    params:
        mem="4000",
        ref=config['reference'],
        prefix="{outpath}/{sample}/02_map/bqsr/{sample}"
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar /home/junhui.li11-umw/anaconda3/envs/snakemake/share/gatk4-4.1.8.1-0/gatk-package-4.1.8.1-local.jar ApplyBQSR -I {input.i1} -O {output.o1} -R {params.ref} --bqsr-recal-file {input.i2} > {log} 2>&1
        samtools flagstat {output.o1} > {params.prefix}.txt
        """