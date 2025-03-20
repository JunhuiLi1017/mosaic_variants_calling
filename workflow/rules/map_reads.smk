rule map_raw:
    input:
        get_clean_fastq
    output:
        o1=temp("{outpath}/02_map/bwa/{sample}/{sample}.raw.bam")
    log:
        l1="{outpath}/02_map/logs/{sample}.bwa.raw.log"
    threads:
        8
    params:
        rg=get_reads_group,
        ref=config['reference']
    resources:
        mem_mb=2000
    shell:
        """
        bwa mem -t {threads} -M {params.rg} {params.ref} {input[0]} {input[1]} | samtools view -@ {threads} -b -o {output.o1} > {log.l1} 2>&1
        """

rule map_raw_sort:
    message:
        "use sambamba instead of samtools sort due to the large mem required by samtools"
    input:
        "{outpath}/02_map/bwa/{sample}/{sample}.raw.bam"
    output:
        sort_bam="{outpath}/02_map/bwa/{sample}/{sample}.sort.bam",
        sort_stat="{outpath}/02_map/bwa/{sample}/{sample}.sort.stat"
    log:
        "{outpath}/02_map/logs/{sample}.raw.sort.log"
    params:
        tmpdir="{outpath}/02_map/bwa/{sample}/tmpdir_{sample}"
    threads:
        8
    resources:
        mem_mb=2000
    shell:
        """
        mkdir -p {params.tmpdir} && \
        sambamba sort \
            -t {threads} \
            -m 2G \
            -o {output.sort_bam} \
            --tmpdir {params.tmpdir} \
            {input} > {log} 2>&1
        samtools stats {output.sort_bam} > {output.sort_stat}
        samtools index {output.sort_bam}
        """

rule remove_dup:
    input:
        "{outpath}/02_map/bwa/{sample}/{sample}.sort.bam"
    output:
        o1=temp("{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.bam"),
        o2="{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.matrix"
    log:
        "{outpath}/02_map/logs/{sample}.sort.rmdup.log"
    threads:
        8
    resources:
        mem_mb=2000
    params:
        dedup="--REMOVE_DUPLICATES true" if config["pcr_based"] == True else "--REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500",
        tmpdir="{outpath}/02_map/dup/{sample}/tmpdir_{sample}"
    shell:
        """
        mkdir -p {params.tmpdir} && picard -Xms10g MarkDuplicates --TMP_DIR {params.tmpdir} --INPUT {input} --METRICS_FILE {output.o2} --OUTPUT {output.o1} {params.dedup} && samtools index {output.o1} > {log} 2>&1
        """

rule BaseRecalibrator:
    input:
        "{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.bam"
    output:
        o1="{outpath}/02_map/bqsr/{sample}/{sample}.recal_data.table"
    log:
        "{outpath}/02_map/logs/{sample}.BQSR.log"
    threads:
        8
    params:
        gatk=config['gatk_current_using'],
        ref=config['reference'],
        dpsnp138=config['dpsnp138'],
        known_indels=config['known_indels'],
        Mills_and_1000G=config['Mills_and_1000G']
    resources:
        mem_mb=2000
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar {params.gatk} BaseRecalibrator -I {input} -O {output.o1} -R {params.ref} --known-sites {params.dpsnp138} --known-sites {params.known_indels} --known-sites {params.Mills_and_1000G} > {log} 2>&1
        """

rule ApplyBQSR:
    input:
        i1="{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.bam",
        i2="{outpath}/02_map/bqsr/{sample}/{sample}.recal_data.table"
    output:
        o1="{outpath}/02_map/bqsr/{sample}/{sample}.sort.rmdup.bqsr.bam",
        o2="{outpath}/02_map/bqsr/{sample}/{sample}.sort.rmdup.bqsr.bai"
    log:
        "{outpath}/02_map/logs/{sample}.applyBQSR.log"
    threads:
        8
    params:
        ref=config['reference'],
        gatk=config['gatk_current_using'],
        prefix="{outpath}/02_map/bqsr/{sample}/{sample}"
    resources:
        mem_mb=2000
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar {params.gatk} ApplyBQSR -I {input.i1} -O {output.o1} -R {params.ref} --bqsr-recal-file {input.i2} > {log} 2>&1
        samtools flagstat {output.o1} > {params.prefix}.txt
        """
