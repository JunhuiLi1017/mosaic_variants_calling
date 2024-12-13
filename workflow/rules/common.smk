import pandas as pd
import glob
import os

from snakemake.utils import validate

configfile: "../config/config.yaml"

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
            f"{outpath}/{wildcards.sample}/01_multiqc/fastp/{wildcards.sample}.R1.fastq.gz",
            f"{outpath}/{wildcards.sample}/01_multiqc/fastp/{wildcards.sample}.R2.fastq.gz"
        ]
    else:
        return [f"{outpath}/{wildcards.sample}/01_multiqc/fastp/{wildcards.sample}.R1.fastq.gz"]
#note: if remove f in f"", an error Error:
#  KeyError: 'Sample {wildcards.sample} not found in units DataFrame.'
#Wildcards:
#  sample={wildcards.sample}
# in get_fastq(wildcards)

def is_single_end(sample):
    """Return True if sample is single end."""
    return pd.isnull(units.loc[(sample), "fq2"])

def get_multiqc_input(wildcards):
    if paired_end:
        GROUP = ["R1", "R2"]
    else:
        GROUP = ["R1"]
    fastqc_files = expand(
        [
            "{outpath}/{sample}/01_multiqc/fastqc/{sample}.{group}_fastqc.html",
            "{outpath}/{sample}/01_multiqc/fastqc/{sample}.{group}_fastqc.zip"
        ],
        sample=[u.sample for u in units.itertuples()],
        group=GROUP
    )

    bqsr_stat_files = expand(
        "{outpath}/{sample}/02_Map/bqsr_stat/{sample}.sort.rmdup.bqsr.stat",
        sample=[u.sample for u in units.itertuples()]
    )
    return fastqc_files + bqsr_stat_files

def get_bqsr_bam(wildcards):
    return expand(["result/02_Map/bqsr/{u.sample}.sort.rmdup.bqsr.bam"], u = units.itertuples())

def get_bqsr_bam_file(wildcards):
    """get fastq files of given sample"""
    bam = units.loc[(wildcards.sample), ["bqsr_bam"]].dropna()
    return {"bqsr_bam": bam.bqsr_bam}

def get_reads_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample), "platform"],
    )

def get_rawbam_summary(wildcards):
    return expand(
            ["{outpath}/{sample}/02_Map/target/{u.sample}.targe.stat.txt"],
             u=units.itertuples()
        )

def get_plot(wildcards):
    return expand(
            ["{outpath}/{sample}/02_Map/target/{u.sample}.bqsr.target.coverage.hist_depth_cumcov.png"],
             u=units.itertuples()
        )