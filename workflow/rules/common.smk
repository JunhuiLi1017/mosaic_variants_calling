import pandas as pd
import glob
import os
import yaml
from snakemake.utils import validate

with open(config["resource"]) as f:
    resource = yaml.safe_load(f)
configfile: "config/config.yaml"
validate(config, schema = "../schemas/config.schema.yaml")
units = pd.read_table(config["units"], dtype = str).set_index(
    ["sample", "library", "flowlane"], drop = False
)
validate(units, schema = "../schemas/units.schema.yaml")

sample_num=units['sample'].unique().tolist()
sample=units['sample']
paired_end = config["paired"]
outpath = config['outpath']
intervals_dir=config["interval"]
hg38_sub_gnomad211_exome_genome_dir=config['hg38_sub_gnomad211_exome_genome_dir']
ref_version=config["ref_version"]
cov=config['depth']

# Get the list of all interval files
interval_files = glob.glob(os.path.join(intervals_dir, "*.intervals.list"))
chromosomes = [os.path.basename(f).replace(".intervals.list", "") for f in interval_files]
standard_order = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
    "chr20", "chr21", "chr22", "chrX", "chrY"
]

chromosomes = sorted(chromosomes, key=lambda x: standard_order.index(x) if x in standard_order else len(standard_order))

def get_fastq(wildcards):
    """Get fastq files for given sample, library, and flowlane."""
    try:
        fastqs = units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), ["fq1", "fq2"]].dropna()
        if len(fastqs) == 2:
            return {"r1": fastqs.fq1, "r2": fastqs.fq2}
        return {"r1": fastqs.fq1}
    except KeyError:
        raise KeyError(f"No entry in units for sample={wildcards.sample}, library={wildcards.library}, flowlane={wildcards.flowlane}")

def get_fastp_fastq(wildcards):
    if paired_end:
        return [
            f"{outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R1.fastq.gz",
            f"{outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R2.fastq.gz"
        ]
    else:
        return [f"{outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R1.fastq.gz"]

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
    return expand(
            [
                "{outpath}/01_multiqc/fastqc/{u.sample}_{u.library}_{u.flowlane}.{group}_fastqc.html",
                "{outpath}/01_multiqc/fastqc/{u.sample}_{u.library}_{u.flowlane}.{group}_fastqc.zip"
            ],
             outpath=outpath,
             u=units.itertuples(),
             group = GROUP
        )

def get_raw_bam(wildcards):
    sample_units = units.loc[units['sample'] == wildcards.sample]
    bam_files = [
        f"{outpath}/02_Map/bwa/sort/{sample}/{row['sample']}_{row['library']}_{row['flowlane']}.sort.bam"
        for _, row in sample_units.iterrows()]
    return bam_files

def get_reads_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{library}_{flowlane}\tSM:{sample}\tPL:{platform}:LB:{library}'".format(
        sample=wildcards.sample,
        #library=units.loc[(wildcards.sample), "library"],
        #flowlane=units.loc[(wildcards.sample, wildcards.library), "flowlane"],
        library=wildcards.library,
        flowlane=wildcards.flowlane,
        platform=units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "platform"]
    )
