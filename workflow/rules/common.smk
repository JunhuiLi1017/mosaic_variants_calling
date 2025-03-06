import pandas as pd
import glob
import os

from snakemake.utils import validate

configfile: "config/config.yaml"

validate(config, schema = "../schemas/config.schema.yaml")

units = pd.read_table(config["units"], dtype = str).set_index(
    ["sample"], drop = False
)
validate(units, schema = "../schemas/units.schema.yaml")

sample=units['sample']

paired_end = config["paired"]

outpath = config['outpath']

trim_value=config['trim_value']

intervals_dir=config["interval"]

ref_version=config["ref_version"]

# Get the list of all interval files
interval_files = glob.glob(os.path.join(intervals_dir, "*.intervals.list"))

# Extract chromosome names by stripping the file path and extension
chromosomes = [os.path.basename(f).replace(".intervals.list", "") for f in interval_files]
cov=config['depth']

def get_fastq(wildcards):
    """get fastq files of given sample"""
    fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_clean_fastq(wildcards):
    if paired_end:
        return [
            f"{outpath}/01_multiqc/fastp/{wildcards.sample}/{wildcards.sample}.R1.fastq.gz",
            f"{outpath}/01_multiqc/fastp/{wildcards.sample}/{wildcards.sample}.R2.fastq.gz"
        ]
    else:
        return [f"{outpath}/01_multiqc/fastp/{wildcards.sample}/{wildcards.sample}.R1.fastq.gz"]
#note: if remove f in f"", an error Error:
#  KeyError: 'Sample {wildcards.sample} not found in units DataFrame.'
#Wildcards:
#  sample={wildcards.sample}
# in get_fastq(wildcards)

def get_multiqc_input(wildcards):
    if paired_end:
        GROUP = ["R1", "R2"]
    else:
        GROUP = ["R1"]
    fastqc_files = expand(
        [
            "{outpath}/01_multiqc/fastqc/{sample}.{group}_fastqc.html",
            "{outpath}/01_multiqc/fastqc/{sample}.{group}_fastqc.zip"
        ],
        sample=[u.sample for u in units.itertuples()],
        group=GROUP
    )

    #bqsr_stat_files = expand(
    #    "{outpath}/02_Map/bqsr_stat/{sample}/{sample}.sort.rmdup.bqsr.stat",
    #    sample=[u.sample for u in units.itertuples()]
    #)
    #return fastqc_files + bqsr_stat_files
    return fastqc_files

def get_reads_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tPL:{platform}\tSM:{sample}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample), "platform"],
    )
