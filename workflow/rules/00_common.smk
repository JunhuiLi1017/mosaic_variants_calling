import pandas as pd
import glob
import os
import yaml
from snakemake.utils import validate

configfile: "config/config.yaml"
validate(config, schema = "../schemas/config.schema.yaml")

with open(config["resource"]) as f:
    resource = yaml.safe_load(f)
validate(resource, schema = "../schemas/resource.schema.yaml")

units = pd.read_table(config["units"], dtype = str).set_index(
    ["sample", "library", "flowlane"], drop = False
)

units['trim_front1'] = pd.to_numeric(units['trim_front1'], errors='coerce').fillna(0).astype(int)
units['trim_front2'] = pd.to_numeric(units['trim_front2'], errors='coerce').fillna(0).astype(int)
units['trim_tail1'] = pd.to_numeric(units['trim_tail1'], errors='coerce').fillna(0).astype(int)
units['trim_tail2'] = pd.to_numeric(units['trim_tail2'], errors='coerce').fillna(0).astype(int)
units['pcr_based'] = units['pcr_based'].map({'Yes': True, 'No': False, 'TRUE': True, 'FALSE': False, '1': True, '0': False, 'True': True, 'False': False}).astype(bool)

validate(units, schema = "../schemas/units.schema.yaml")


sample_num=units['sample'].unique().tolist()
Sample=units['sample']
Outpath = config['outpath']
intervals_dir=config["interval"]
Ref_version=config["ref_version"]
Callers=config["callers"]
Coverage=config['coverage']
Model=config['mosaic_pred_model']
hg38_sub_gnomad211_exome_genome_dir=config['hg38_sub_gnomad211_exome_genome_dir']
# Get the list of all interval files
interval_files = glob.glob(os.path.join(intervals_dir, "*.intervals.list"))
CHROMOSOMES = [os.path.basename(f).replace(".intervals.list", "") for f in interval_files]
standard_order = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
    "chr20", "chr21", "chr22", "chrX", "chrY"
]

CHROMOSOMES = sorted(CHROMOSOMES, key=lambda x: standard_order.index(x) if x in standard_order else len(standard_order))

def is_pair_end(wildcards):
    try:
        fastqs = units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), ["fq1", "fq2"]].dropna()
        if len(fastqs) == 2:
            return True
        return False
    except KeyError:
        raise KeyError(f"No entry in units for sample={wildcards.sample}, library={wildcards.library}, flowlane={wildcards.flowlane}")

def get_fastq(wildcards):
    """Get fastq files for given sample, library, and flowlane."""
    fastqs = units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), ["fq1", "fq2"]].dropna()
    if is_pair_end:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_fastqc_fastp(wildcards):
    if is_pair_end:
        return [
            f"{Outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R1.fastq.gz",
            f"{Outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R2.fastq.gz"
        ]
    else:
        return [f"{Outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R1.fastq.gz"]

def get_multiqc_fastp(wildcards):
    if is_pair_end:
        GROUP = ["R1", "R2"]
    else:
        GROUP = ["R1"]
    return expand(
            [
                f"{Outpath}/01_multiqc/fastqc/{{sample}}_{{library}}_{{flowlane}}.{{group}}_fastqc.html",
                f"{Outpath}/01_multiqc/fastqc/{{sample}}_{{library}}_{{flowlane}}.{{group}}_fastqc.zip"
            ],
            sample=[u.sample for u in units.itertuples()],
            library=[u.library for u in units.itertuples()],
            flowlane=[u.flowlane for u in units.itertuples()],
            group = GROUP
        )

def get_vcf_inputs(wildcards):
    callers = config["callers"]
    input_list = []
    for caller in callers:
        input_list.append(f"{Outpath}/03_variants/{caller}/04_pass/{wildcards.sample}.{Ref_version}.pass.vcf.gz")
    return input_list

def get_raw_bam(wildcards):
    sample_units = units.loc[units['sample'] == wildcards.sample]
    bam_files = [
        f"{Outpath}/02_Map/bwa/sort/{sample}/{row['sample']}_{row['library']}_{row['flowlane']}.sort.bam"
        for _, row in sample_units.iterrows()]
    return bam_files

def get_reads_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{library}_{flowlane}\tSM:{sample}\tPL:{platform}\tLB:{library}'".format(
        sample=wildcards.sample,
        library=wildcards.library,
        flowlane=wildcards.flowlane,
        platform=units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "platform"]
    )