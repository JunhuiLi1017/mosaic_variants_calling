rule fastp:
    input:
        unpack(get_fastq)
    output:
        temp(["{outpath}/01_multiqc/fastp/{sample}/{sample}.R1.fastq.gz", "{outpath}/01_multiqc/fastp/{sample}/{sample}.R2.fastq.gz"]) if paired_end else temp("{outpath}/01_multiqc/fastp/{sample}/{sample}.R1.fastq.gz"),
        j="{outpath}/01_multiqc/fastp/{sample}/{sample}.json", 
        h="{outpath}/01_multiqc/fastp/{sample}/{sample}.html"
    log:
        "{outpath}/01_multiqc/logs/{sample}.fastp.log"
    params:
        trim_expr= "" if config['trim'] == False else f"-f {trim_value} -F {trim_value}",
        in_r2=lambda wildcards, input: f"-I {input.r2}" if paired_end else "",
        out_r2=lambda wildcards, output: f"-O {output[1]}" if paired_end else ""
    threads: 4  # Increased to leverage fastp multi-threading
    resources:
        mem_mb=1500
    shell:
        """
        fastp {params.trim_expr} \
            -i {input.r1} \
            {params.in_r2} \
            -o {output[0]} \
            {params.out_r2} \
            -j {output.j} -h {output.h} \
            --thread {threads} \
            > {log} 2>&1
        """

rule fastqc:
    input:
        get_clean_fastq
    output:
        ["{outpath}/01_multiqc/fastqc/{sample}.R1_fastqc.html", "{outpath}/01_multiqc/fastqc/{sample}.R2_fastqc.html"] if paired_end else "{outpath}/01_multiqc/fastqc/{sample}.R1_fastqc.html",
        ["{outpath}/01_multiqc/fastqc/{sample}.R1_fastqc.zip", "{outpath}/01_multiqc/fastqc/{sample}.R2_fastqc.zip"] if paired_end else "{outpath}/01_multiqc/fastqc/{sample}.R1_fastqc.zip"
    conda:
        "../envs/multiqc_env.yaml"
    log:
        "{outpath}/01_multiqc/logs/{sample}.fastqc.log"
    params:
        out_fastqc="{outpath}/01_multiqc/fastqc/{sample}"
    resources:
        mem_mb=1000
    shell:
        "fastqc -o {params.out_fastqc} {input} > {log} 2>&1"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "{outpath}/01_multiqc/multiqc_report.html"
    params:
        outdir="{outpath}/01_multiqc",
        indir="{outpath}/01_multiqc/fastqc"
    conda:
        "../envs/multiqc_env.yaml"
    log:
        "{outpath}/01_multiqc/logs/multiqc.log"
    resources:
        mem_mb=1000
    shell:
        """
        multiqc -o {params.outdir} {params.indir} --force
        """