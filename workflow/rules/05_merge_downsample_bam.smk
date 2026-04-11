
"""
Workflow for merging CX/STR BAM files, downsampling, and read group replacement.
This workflow processes BAM files that need to be merged, downsampled, and have read groups replaced.

Configuration:
- Provide a TSV file in config['merge_downsample_units'] with columns:
  - sample: sample identifier (e.g., 'sample1', 'NC1', etc.)
  - cx_bam: path to CX BAM file
  - str_bam: path to STR BAM file
  - downsample_seed: (optional) 6-digit random seed for downsampling
    If not provided, will be generated deterministically from sample name.
"""

import os
import pandas as pd
import re
import hashlib

# Read merge/downsample units from config
merge_downsample_units_file = config.get('merge_downsample_units', 'config/merge_downsample_units.tsv')
if os.path.exists(merge_downsample_units_file):
    merge_downsample_units = pd.read_table(merge_downsample_units_file, dtype=str)
    # Ensure downsample_seed column exists, generate if missing
    if 'downsample_seed' not in merge_downsample_units.columns:
        merge_downsample_units['downsample_seed'] = None
    # Generate downsample_seed for rows where it's missing
    for idx, row in merge_downsample_units.iterrows():
        if pd.isna(row.get('downsample_seed')) or str(row.get('downsample_seed', '')).strip() == '':
            seed_hash = int(hashlib.md5(row['sample'].encode()).hexdigest()[:6], 16) % 1000000
            merge_downsample_units.at[idx, 'downsample_seed'] = str(seed_hash).zfill(6)
else:
    merge_downsample_units = pd.DataFrame(columns=['sample', 'cx_bam', 'str_bam', 'downsample_seed'])

def extract_prefix_and_seed(sample_name, cx_bam_path):
    """
    Extract prefix and seed from sample name or BAM filename.
    Tries sample name first (e.g., 'NC1-sampleseed_123'), then BAM filename.
    """
    # Try to extract from sample name
    match = re.search(r'([^-]+)-sampleseed_(\d+)', sample_name)
    if match:
        return match.group(1), match.group(2)
    
    # Try to extract from CX BAM filename
    bam_basename = os.path.basename(cx_bam_path)
    match = re.search(r'([^-]+)-CX.*sampleseed_(\d+)', bam_basename)
    if match:
        return match.group(1), match.group(2)
    
    # Try simpler pattern: prefix at start, seed anywhere
    match = re.search(r'^([^-]+)', sample_name)
    prefix = match.group(1) if match else sample_name.split('-')[0]
    
    seed_match = re.search(r'sampleseed_(\d+)', bam_basename)
    seed = seed_match.group(1) if seed_match else '0'
    
    return prefix, seed

def get_sample_info(sample_name):
    """Get all info for a sample from the units DataFrame."""
    row = merge_downsample_units[merge_downsample_units['sample'] == sample_name]
    if row.empty:
        raise ValueError(f"Sample {sample_name} not found in merge_downsample_units")
    return row.iloc[0]

def generate_downsample_seed(sample_name):
    """Generate deterministic downsample seed from sample name."""
    seed_hash = int(hashlib.md5(sample_name.encode()).hexdigest()[:6], 16) % 1000000
    return str(seed_hash).zfill(6)

def get_cx_bam(wildcards):
    """Get CX BAM file for given sample."""
    row = get_sample_info(wildcards.sample)
    return row['cx_bam']

def get_str_bam(wildcards):
    """Get STR BAM file for given sample."""
    row = get_sample_info(wildcards.sample)
    return row['str_bam']

def get_prefix(wildcards):
    """Get prefix for given sample."""
    row = get_sample_info(wildcards.sample)
    prefix, _ = extract_prefix_and_seed(wildcards.sample, row['cx_bam'])
    return prefix

def get_seed(wildcards):
    """Get seed for given sample."""
    row = get_sample_info(wildcards.sample)
    _, seed = extract_prefix_and_seed(wildcards.sample, row['cx_bam'])
    return seed

def get_downsample_seed(wildcards):
    """Get downsample seed for given sample and downsample_seed wildcard."""
    row = get_sample_info(wildcards.sample)
    seed = str(row['downsample_seed']).strip()
    # Verify it matches the wildcard
    if seed != str(wildcards.downsample_seed):
        raise ValueError(f"Downsample seed mismatch: config has {seed}, wildcard has {wildcards.downsample_seed}")
    return seed

rule merge_cx_str_bam:
    """
    Merge CX and STR BAM files for each sample.
    """
    input:
        cx_bam=get_cx_bam,
        str_bam=get_str_bam
    output:
        merged_bam="{outpath}/05_merge_downsample/01_merged/{sample}.merged.sort.rmdup.bqsr.bam"
    params:
        prefix=get_prefix,
        seed=get_seed
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    log:
        "{outpath}/05_merge_downsample/logs/{sample}.merge.log"
    container:
        container_image["samtools_1.20"]
    shell:
        """
        samtools merge -@ {threads} {output.merged_bam} {input.cx_bam} {input.str_bam} > {log} 2>&1
        """

rule downsample_bam:
    """
    Downsample merged BAM files to 50% using a specified seed.
    The downsample_seed can be provided in the config or will be generated deterministically.
    """
    input:
        "{outpath}/05_merge_downsample/01_merged/{sample}.merged.sort.rmdup.bqsr.bam"
    output:
        downsampled_bam=temp("{outpath}/05_merge_downsample/02_downsampled/{sample}.{downsample_seed}.50pct.tmp.bam"),
        sorted_bam="{outpath}/05_merge_downsample/02_downsampled/{sample}.{downsample_seed}.50pct.sorted.bam",
        bai="{outpath}/05_merge_downsample/02_downsampled/{sample}.{downsample_seed}.50pct.sorted.bam.bai"
    params:
        downsample_seed=get_downsample_seed,
        fraction="0.5"
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    log:
        "{outpath}/05_merge_downsample/logs/{sample}.{downsample_seed}.50pct.log"
    container:
        container_image["samtools_1.20"]
    shell:
        """
        samtools view -@ {threads} -b -s {params.downsample_seed}{params.fraction} {input} > {output.downsampled_bam} 2>> {log}
        samtools sort -@ {threads} -o {output.sorted_bam} {output.downsampled_bam} >> {log} 2>&1
        samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1
        rm -f {output.downsampled_bam}
        """

rule replace_read_groups:
    """
    Replace read groups in downsampled BAM files with a single new sample name.
    Sample ID format: {prefix}.{downsample_seed}
    """
    input:
        bam="{outpath}/05_merge_downsample/02_downsampled/{sample}.{downsample_seed}.50pct.sorted.bam",
        bai="{outpath}/05_merge_downsample/02_downsampled/{sample}.{downsample_seed}.50pct.sorted.bam.bai"
    output:
        bam="{outpath}/05_merge_downsample/03_rg_replaced/{sample}.{downsample_seed}.50pct.sorted.bam.AddOrReplaceReadGroups.bam",
        bai="{outpath}/05_merge_downsample/03_rg_replaced/{sample}.{downsample_seed}.50pct.sorted.bam.AddOrReplaceReadGroups.bam.bai"
    params:
        prefix=get_prefix,
        downsample_seed=get_downsample_seed,
        sample_id=lambda wildcards: f"{get_prefix(wildcards)}.{get_downsample_seed(wildcards)}",
        read_group_id=lambda wildcards: f"{get_prefix(wildcards)}.{get_downsample_seed(wildcards)}",
        library="A",
        platform="ILLUMINA",
        platform_unit="1",
        sample_name=lambda wildcards: f"{get_prefix(wildcards)}.{get_downsample_seed(wildcards)}",
        picard_jar=config.get('picard_jar', '/home/junhui.li11-umw/anaconda3/envs/picard3/share/picard-3.0.0-1/picard.jar'),
        java_mem=lambda wildcards, resources: f"-Xmx{resources.mem_mb // 1000}g"
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    log:
        "{outpath}/05_merge_downsample/logs/{sample}.{downsample_seed}.replace_rg.log"
    container:
        container_image.get("picard_3.0.0", "picard:3.0.0")
    shell:
        """
        java {params.java_mem} -XX:ParallelGCThreads={threads} \
            -jar {params.picard_jar} AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -ID {params.read_group_id} \
            -LB {params.library} \
            -PL {params.platform} \
            -PU {params.platform_unit} \
            -SM {params.sample_name} > {log} 2>&1
        samtools index {output.bam}
        mv {output.bam}.bai {output.bai}
        """

rule create_symlinks:
    """
    Create symlinks for processed BAM files with simplified naming.
    Sample ID format: {prefix}.{downsample_seed}
    """
    input:
        bam="{outpath}/05_merge_downsample/03_rg_replaced/{sample}.{downsample_seed}.50pct.sorted.bam.AddOrReplaceReadGroups.bam",
        bai="{outpath}/05_merge_downsample/03_rg_replaced/{sample}.{downsample_seed}.50pct.sorted.bam.AddOrReplaceReadGroups.bam.bai"
    output:
        bam_link="{outpath}/05_merge_downsample/04_symlinks/{sample_id}.bam",
        bai_link="{outpath}/05_merge_downsample/04_symlinks/{sample_id}.bai"
    params:
        sample_id=lambda wildcards: f"{get_prefix(wildcards)}.{get_downsample_seed(wildcards)}"
    shell:
        """
        ln -sf $(readlink -f {input.bam}) {output.bam_link}
        ln -sf $(readlink -f {input.bai}) {output.bai_link}
        """

def get_all_sample_ids():
    """Get all sample IDs from merge_downsample_units."""
    sample_ids = []
    for _, row in merge_downsample_units.iterrows():
        sample_name = row['sample']
        prefix, _ = extract_prefix_and_seed(sample_name, row['cx_bam'])
        if pd.notna(row['downsample_seed']) and str(row['downsample_seed']).strip():
            downsample_seed = str(row['downsample_seed']).strip()
        else:
            downsample_seed = generate_downsample_seed(sample_name)
        sample_ids.append(f"{prefix}.{downsample_seed}")
    return sample_ids

rule generate_units_tsv:
    """
    Generate units.tsv file for all processed samples.
    Uses the merge_downsample_units config to determine all sample IDs.
    """
    input:
        bam_files=expand(
            "{outpath}/05_merge_downsample/04_symlinks/{sample_id}.bam",
            sample_id=get_all_sample_ids()
        )
    output:
        units_tsv="{outpath}/05_merge_downsample/units_downsampled.tsv"
    params:
        depth=config.get('downsample_depth', '50'),
        sex=config.get('downsample_sex', 'M'),
        platform="ILLUMINA",
        library="A",
        flowlane="2201"
    shell:
        """
        echo -e "sample\\tlibrary\\tflowlane\\tplatform\\tfq1\\tfq2\\tdepth\\tsex\\ttrim_front1\\ttrim_front2\\ttrim_tail1\\ttrim_tail2\\tpcr_based" > {output.units_tsv}
        for bam in {input.bam_files}; do
            if [ -f "$bam" ]; then
                sample=$(basename "$bam" .bam)
                echo -e "$sample\\t{params.library}\\t{params.flowlane}\\t{params.platform}\\t$sample.fq1.gz\\t$sample.fq2.gz\\t{params.depth}\\t{params.sex}\\t0\\t0\\t0\\t0\\tFalse" >> {output.units_tsv}
            fi
        done
        """

