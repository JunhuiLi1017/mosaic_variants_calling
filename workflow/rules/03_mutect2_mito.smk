rule subset_bam_to_chrm:
    input:
        bam = "{outpath}/02_map/02_sort/{sample}.{library}.{flowlane}.sort.bam",      # or .cram
        bai = "{outpath}/02_map/02_sort/{sample}.{library}.{flowlane}.sort.bam.bai"   # or .cram.crai
    output:
        bam = "{outpath}/02_map/00_chrM/01_subset_to_chrm/{sample}.{library}.{flowlane}.sort.chrm.bam"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.{library}.{flowlane}.subset_to_chrm.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    params:
        tmpdir="{outpath}/02_map/00_chrM/01_subset_to_chrm/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000) // threads,
        ref = config['reference_chrM']
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" PrintReads \
            -L chrM \
            -I {input.bam} \
            -O {output.bam} \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            --tmp-dir {params.tmpdir} \
            > {log} 2>&1
        """

rule merge_chrm_bams:
    input:
        lambda wildcards: [
            f"{wildcards.outpath}/02_map/00_chrM/01_subset_to_chrm/{wildcards.sample}.{u.library}.{u.flowlane}.sort.chrm.bam"
            for u in units[units["sample"] == wildcards.sample].itertuples()
        ]
    output:
        bam = "{outpath}/02_map/00_chrM/02_merge_chrm_bams/{sample}.merged.chrm.bam",
        bai = "{outpath}/02_map/00_chrM/02_merge_chrm_bams/{sample}.merged.chrm.bam.bai"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.merge_chrm_bams.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    container:
        container_image["samtools_1.20"]
    shell:
        """
        samtools merge -@ {threads} -O BAM {output.bam} {input} > {log} 2>&1
        samtools index -@ {threads} {output.bam} > {log} 2>&1
        """

rule revert_sam:
    input:
        bam = "{outpath}/02_map/00_chrM/02_merge_chrm_bams/{sample}.merged.chrm.bam"   # your aligned input BAM
    output:
        bam = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.bam"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.revert_sam.log"
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    params:
        tmpdir="{outpath}/02_map/00_chrM/03_unmapped_chrm/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:    
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" RevertSam \
            -I {input.bam} \
            -O {output.bam} \
            --OUTPUT_BY_READGROUP false \
            --VALIDATION_STRINGENCY LENIENT \
            --ATTRIBUTE_TO_CLEAR FT \
            --ATTRIBUTE_TO_CLEAR CO \
            --SORT_ORDER queryname \
            --RESTORE_ORIGINAL_QUALITIES false \
            > {log} 2>&1
        """

rule bam_to_fq:
    input:
        bam = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.bam"
    output:
        fq1 = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.1.fq",
        fq2 = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.2.fq"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.bam_to_fq.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    container:
        container_image["bedtools_2.27.1"]
    shell:
        """
        bedtools bamtofastq -i {input.bam} -fq {output.fq1} -fq2 {output.fq2} > {log} 2>&1
        """

rule align_reads:
    input:
        fq1 = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.1.fq",
        fq2 = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.2.fq"
    output:
        bam = "{outpath}/02_map/00_chrM/04_align/{sample}.merged.chrm.unmapped.aligned.bam"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.align_reads.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        ref = config['reference_chrM']
    shell:
        """
        bwa mem -K 100000000 -v 3 -t {threads} -Y {params.ref} \
            {input.fq1} {input.fq2} \
            | samtools view -@ {threads} -b -o {output.bam} - 2> {log}
        """

rule align_reads_shifted:
    input:
        fq1 = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.1.fq",
        fq2 = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.2.fq"
    output:
        bam = "{outpath}/02_map/00_chrM/04_align_shifted/{sample}.merged.chrm.unmapped.aligned.shifted.bam"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.align_reads_shifted.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        ref = config['reference_chrM_shifted']
    shell:
        """
        bwa mem -K 100000000 -v 3 -t {threads} -Y {params.ref} \
            {input.fq1} {input.fq2} \
            | samtools view -@ {threads} -b -o {output.bam} - 2> {log}
        """

rule merge_bams_aligned:
    input:
        bams = "{outpath}/02_map/00_chrM/04_align/{sample}.merged.chrm.unmapped.aligned.bam",
        unmap_bam = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.bam"
    output:
        bam = "{outpath}/02_map/00_chrM/05_merge_unmap_map/{sample}.merge_unmap_map.bam"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.merge_unmap_map.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        ref = config['reference_chrM'],
        tmpdir="{outpath}/02_map/00_chrM/05_merge_unmap_map/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" MergeBamAlignment \
            --ALIGNED_BAM {input.bams} \
            --UNMAPPED_BAM {input.unmap_bam} \
            --OUTPUT {output.bam} \
            --REFERENCE_SEQUENCE {params.ref} \
            --ADD_MATE_CIGAR true \
            --CLIP_ADAPTERS false \
            --TMP_DIR {params.tmpdir} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ATTRIBUTES_TO_REMOVE NM \
            --ATTRIBUTES_TO_REMOVE MD \
            --SORT_ORDER "unsorted" \
            --IS_BISULFITE_SEQUENCE false \
            --ALIGNED_READS_ONLY false \
            --MAX_RECORDS_IN_RAM 2000000 \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --PROGRAM_RECORD_ID "bwamem" \
            --PROGRAM_GROUP_VERSION "0.7.17" \
            --PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -v 3 -t {threads} -Y {params.ref}" \
            > {log} 2>&1
        """

rule merge_bams_aligned_shifted:
    input:
        bams = "{outpath}/02_map/00_chrM/04_align_shifted/{sample}.merged.chrm.unmapped.aligned.shifted.bam",
        unmap_bam = "{outpath}/02_map/00_chrM/03_unmapped_chrm/{sample}.merged.chrm.unmapped.bam"
    output:
        bam = "{outpath}/02_map/00_chrM/05_merge_unmap_map_shifted/{sample}.merge_unmap_map_shifted.bam"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.merge_unmap_map_shifted.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        ref = config['reference_chrM_shifted'],
        tmpdir="{outpath}/02_map/00_chrM/05_merge_unmap_map_shifted/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" MergeBamAlignment \
            --ALIGNED_BAM {input.bams} \
            --UNMAPPED_BAM {input.unmap_bam} \
            --OUTPUT {output.bam} \
            --REFERENCE_SEQUENCE {params.ref} \
            --ADD_MATE_CIGAR true \
            --CLIP_ADAPTERS false \
            --TMP_DIR {params.tmpdir} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ATTRIBUTES_TO_REMOVE NM \
            --ATTRIBUTES_TO_REMOVE MD \
            --SORT_ORDER "unsorted" \
            --IS_BISULFITE_SEQUENCE false \
            --ALIGNED_READS_ONLY false \
            --MAX_RECORDS_IN_RAM 2000000 \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --PROGRAM_RECORD_ID "bwamem" \
            --PROGRAM_GROUP_VERSION "0.7.17" \
            --PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -v 3 -t {threads} -Y {params.ref}" \
            > {log} 2>&1
        """

rule sort_bam:
    input:
        bam = "{outpath}/02_map/00_chrM/05_merge_unmap_map/{sample}.merge_unmap_map.bam"
    output:
        bam = "{outpath}/02_map/00_chrM/06_sort_bam/{sample}.merge_unmap_map.sorted.bam",
        bai = "{outpath}/02_map/00_chrM/06_sort_bam/{sample}.merge_unmap_map.sorted.bam.bai"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.sort_bam.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    params:
        tmpdir="{outpath}/02_map/00_chrM/06_sort_bam/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["sambamba_1.0.1"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        sambamba sort -t {threads} -m {params.command_mem}M -o {output.bam} --tmpdir {params.tmpdir} {input.bam} > {log} 2>&1
        sambamba index -t {threads} {output.bam} > {log} 2>&1
        """

rule sort_bam_shifted:
    input:
        bam = "{outpath}/02_map/00_chrM/05_merge_unmap_map_shifted/{sample}.merge_unmap_map_shifted.bam"
    output:
        bam = "{outpath}/02_map/00_chrM/06_sort_bam_shifted/{sample}.merge_unmap_map_shifted.sorted.bam",
        bai = "{outpath}/02_map/00_chrM/06_sort_bam_shifted/{sample}.merge_unmap_map_shifted.sorted.bam.bai"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.sort_bam_shifted.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    params:
        tmpdir="{outpath}/02_map/00_chrM/06_sort_bam_shifted/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["sambamba_1.0.1"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        sambamba sort -t {threads} -m {params.command_mem}M -o {output.bam} --tmpdir {params.tmpdir} {input.bam} > {log} 2>&1
        sambamba index -t {threads} {output.bam} > {log} 2>&1
        """

rule mark_duplicates:
    input:
        bam = "{outpath}/02_map/00_chrM/06_sort_bam/{sample}.merge_unmap_map.sorted.bam"
    output:
        bam = "{outpath}/02_map/00_chrM/07_mark_duplicates/{sample}.merge_unmap_map.sorted.rmdup.bam",
        metrics = "{outpath}/02_map/00_chrM/07_mark_duplicates/{sample}.merge_unmap_map.sorted.rmdup.metrics.txt"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.markdup_bam.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    params:
        tmpdir="{outpath}/02_map/00_chrM/07_mark_duplicates/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --REMOVE_SEQUENCING_DUPLICATES true \
            --CREATE_INDEX true \
            --TMP_DIR {params.tmpdir} \
            > {log} 2>&1
        """

rule mark_duplicates_shifted:
    input:
        bam = "{outpath}/02_map/00_chrM/06_sort_bam_shifted/{sample}.merge_unmap_map_shifted.sorted.bam"
    output:
        bam = "{outpath}/02_map/00_chrM/07_mark_duplicates_shifted/{sample}.merge_unmap_map_shifted.sorted.rmdup.bam",
        metrics = "{outpath}/02_map/00_chrM/07_mark_duplicates_shifted/{sample}.merge_unmap_map_shifted.sorted.rmdup.metrics.txt"
    log:
        "{outpath}/02_map/00_chrM/logs/{sample}.markdup_bam_shifted.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    params:
        tmpdir="{outpath}/02_map/00_chrM/07_mark_duplicates_shifted/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --REMOVE_SEQUENCING_DUPLICATES true \
            --CREATE_INDEX true \
            --TMP_DIR {params.tmpdir} \
            > {log} 2>&1
        """

rule mutect2_mito:
    input:
        bam = "{outpath}/02_map/00_chrM/07_mark_duplicates/{sample}.merge_unmap_map.sorted.rmdup.bam"
    output:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/01_regular/{sample}.vcf.gz",
        idx = "{outpath}/03_variants/mutect2/00_chrM/01_regular/{sample}.vcf.gz.tbi",
        stats = "{outpath}/03_variants/mutect2/00_chrM/01_regular/{sample}.vcf.gz.stats"
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.mutect2_mito.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads,
        extra_args = "--annotation StrandBiasBySample --mitochondria-mode --max-reads-per-alignment-start 75 --max-mnp-distance 0",
        m2_extra_args = "-L chrM:576-16024",
        ref = config['reference_chrM']
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" Mutect2 \
            -R {params.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            {params.extra_args} \
            {params.m2_extra_args} \
            > {log} 2>&1
        """

rule mutect2_mito_shifted:
    input:
        bam = "{outpath}/02_map/00_chrM/07_mark_duplicates_shifted/{sample}.merge_unmap_map_shifted.sorted.rmdup.bam"
    output:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/02_shifted/{sample}.vcf.gz",
        idx = "{outpath}/03_variants/mutect2/00_chrM/02_shifted/{sample}.vcf.gz.tbi",
        stats = "{outpath}/03_variants/mutect2/00_chrM/02_shifted/{sample}.vcf.gz.stats"
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.mutect2_mito_shifted.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    params:
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads,
        extra_args = "--annotation StrandBiasBySample --mitochondria-mode --max-reads-per-alignment-start 75 --max-mnp-distance 0",
        m2_extra_args = "-L chrM:8025-9144",
        ref = config['reference_chrM_shifted']
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" Mutect2 \
            -R {params.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            {params.extra_args} \
            {params.m2_extra_args} \
            > {log} 2>&1
        """

rule merge_stats:
    input:
        regular_stats = "{outpath}/03_variants/mutect2/00_chrM/01_regular/{sample}.vcf.gz.stats",
        shifted_stats = "{outpath}/03_variants/mutect2/00_chrM/02_shifted/{sample}.vcf.gz.stats"
    output:
        stats = "{outpath}/03_variants/mutect2/00_chrM/04_combine/{sample}.final.vcf.gz.stats"
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.merge_stats.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    params:
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" MergeMutectStats \
            --stats {input.shifted_stats} \
            --stats {input.regular_stats} \
            -O {output.stats} \
            > {log} 2>&1
        """

rule liftover_vcfs:
    input:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/02_shifted/{sample}.vcf.gz"
    output:
        lifted_vcf = "{outpath}/03_variants/mutect2/00_chrM/03_liftover/{sample}.shifted.liftover.vcf.gz",
        rejected_vcf = "{outpath}/03_variants/mutect2/00_chrM/03_liftover/{sample}.shifted.liftover.rejected.vcf.gz"
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.liftover_vcfs.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        ref = config['reference_chrM'],
        chain = config['shift_back_chain'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" LiftoverVcf \
            --INPUT {input.vcf} \
            --OUTPUT {output.lifted_vcf} \
            --REFERENCE_SEQUENCE {params.ref} \
            --CHAIN {params.chain} \
            --REJECT {output.rejected_vcf} \
            > {log} 2>&1
        """

rule combine_vcfs:
    input:
        shifted_vcf = "{outpath}/03_variants/mutect2/00_chrM/03_liftover/{sample}.shifted.liftover.vcf.gz",
        regular_vcf = "{outpath}/03_variants/mutect2/00_chrM/01_regular/{sample}.vcf.gz"
    output:
        final_vcf = "{outpath}/03_variants/mutect2/00_chrM/04_combine/{sample}.final.vcf.gz",
        final_idx = "{outpath}/03_variants/mutect2/00_chrM/04_combine/{sample}.final.vcf.gz.tbi"
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.combine_vcfs.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000) // threads
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" MergeVcfs \
            -I {input.shifted_vcf} \
            -I {input.regular_vcf} \
            -O {output.final_vcf} \
            --CREATE_INDEX true \
            > {log} 2>&1
        """

rule get_contamination:
    input:
        bam = "{outpath}/02_map/00_chrM/07_mark_duplicates/{sample}.merge_unmap_map.sorted.rmdup.bam"
    output:
        major_hg      = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.major_hg.txt",
        major_level   = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.major_level.txt",
        minor_hg      = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.minor_hg.txt",
        minor_level   = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.minor_level.txt"
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.get_contamination.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb = resource['resource']['high']['mem_mb']
    params:
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000) // threads,
        ref = config['reference_chrM'],
        outpath = "{outpath}/03_variants/mutect2/00_chrM/05_filter/haplochecker_out",
        vaf= config['contamination_vaf']
    container:
        container_image["mtdnaserver_1.2"]
    shell:
        """
        java -jar /usr/mtdnaserver/mitolib.jar haplochecker \
            --in {input.bam} \
            --ref {params.ref} \
            --out {params.outpath} \
            --QUAL 20 \
            --MAPQ 30 \
            --VAF {params.vaf} > {log} 2>&1

        awk -F'\t' 'NR==2 {{
            print $3 > "{output.major_hg}"
            print $4 > "{output.major_level}"
            print $7 > "{output.minor_hg}"
            print $8 > "{output.minor_level}"
        }}' {params.outpath}/{wildcards.sample}.merge_unmap_map.sorted.rmdup.contamination.txt > {log} 2>&1
        """

rule filter_mutect_calls:
    input:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/04_combine/{sample}.final.vcf.gz",
        stats = "{outpath}/03_variants/mutect2/00_chrM/04_combine/{sample}.final.vcf.gz.stats",
        minor_level = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.minor_level.txt"
    output:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.filtered.vcf.gz",
        idx = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.filtered.vcf.gz.tbi"
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.filter_mutect_calls.log"
    params:
        ref = config['reference_chrM'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000) // threads,
        min_af = config['filter_mutect_calls_min_af'],
        autosomal_coverage = 593
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" FilterMutectCalls \
            -V {input.vcf} \
            -R {params.ref} \
            -O {output.vcf} \
            --stats {input.stats} \
            --mitochondria-mode \
            --min-allele-fraction {params.min_af} \
            --contamination-estimate $(cat {input.minor_level}) > {log} 2>&1
        """

rule variant_filtration:
    input:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.filtered.vcf.gz",
        idx = "{outpath}/03_variants/mutect2/00_chrM/05_filter/{sample}.filtered.vcf.gz.tbi"
    output:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/06_variant_filtration/{sample}.variant_filtration.vcf.gz",
        idx = "{outpath}/03_variants/mutect2/00_chrM/06_variant_filtration/{sample}.variant_filtration.vcf.gz.tbi"
    params:
        blacklist = config['blacklist_site'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000) // threads
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.variant_filtration.log"
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    container:
        container_image["gatk_4.6.1.0"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" VariantFiltration \
            -V {input.vcf} \
            -O {output.vcf} \
            --mask {params.blacklist} \
            --mask-name "blacklisted_site" > {log} 2>&1
        """

rule variant_pass:
    input:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/06_variant_filtration/{sample}.variant_filtration.vcf.gz",
        idx = "{outpath}/03_variants/mutect2/00_chrM/06_variant_filtration/{sample}.variant_filtration.vcf.gz.tbi"
    output:
        vcf = "{outpath}/03_variants/mutect2/00_chrM/07_variant_pass/{sample}.variant_pass.vcf.gz"
    params:
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000) // threads
    log:
        "{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.variant_pass.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    container:
        container_image["bcftools_1.9"]
    shell:
        """
        bcftools view -f PASS {input.vcf} -Oz -o {output.vcf} 2> {log}
        """

rule annotate_clinvar:
	input:
		vcf = "{outpath}/03_variants/mutect2/00_chrM/07_variant_pass/{sample}.variant_pass.vcf.gz"
	output:
		txt="{outpath}/03_variants/mutect2/00_chrM/08_annovar/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/mutect2/00_chrM/logs/{sample}.{ref_version}.annotate_clinvar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/mutect2/00_chrM/08_annovar/{sample}"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["terra_perl_anno"]
	shell:
		"""
		perl {params.annovar_dir}/table_annovar.pl \
		{input.vcf} \
		{params.annovar_dir}/humandb_{params.ref_version} \
		-buildver {params.ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol ensGene,dbnsfp42a,clinvar_20240917,gnomad41_genome,gnomad41_exome \
		-operation g,f,f,f,f \
		-nastring . \
		-vcfinput
		"""