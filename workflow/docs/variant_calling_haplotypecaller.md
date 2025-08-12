# HaplotypeCaller Workflow Documentation

This document describes the HaplotypeCaller workflow implemented in `rules/variant_calling_haplotypecaller.smk`. The workflow uses GATK's HaplotypeCaller for variant calling.

## Overview

The workflow processes BAM files to identify variants using GATK's HaplotypeCaller, which can detect both SNPs and indels with high accuracy.

## Rules

### 1. haplotypecaller
- **Purpose**: Performs variant calling using HaplotypeCaller
- **Input**: BAM file after BQSR
- **Output**: 
  - VCF file with variants
  - Index file
- **Key Parameters**:
  - Reference genome
  - dbSNP database
  - PCR indel model settings
  - Emit reference confidence (GVCF)
- **Tools**: GATK 4.6.1.0

### 2. mergevcfs
- **Purpose**: Merges chromosome-specific VCF files
- **Input**: 
  - Individual chromosome VCF files
  - Index files
- **Output**: 
  - Merged VCF file
  - Index file
- **Tools**: GATK MergeVcfs

### 3. haplotypecaller_anno
- **Purpose**: Performs initial annotation of variants
- **Input**: 
  - VCF file
  - Index file
- **Output**: 
  - Input annotation file
  - VCF annotation file
- **Tools**: ANNOVAR

## Dependencies

The workflow requires several external tools:
- GATK 4.6.1.0
- ANNOVAR
- Reference genome
- dbSNP database

## Configuration

The workflow uses several configuration parameters:
- Reference genome
- dbSNP database
- PCR-based vs non-PCR-based settings
- Resource allocation settings

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
└── 03_variants/
    └── 02_haplotypecaller/
        ├── sub_vcf/
        │   └── {sample}/
        │       ├── {sample}.{chr}.vcf.gz
        │       └── {sample}.{chr}.vcf.gz.tbi
        ├── {sample}/
        │   └── {sample}.vcf.gz
        └── sub_anno/
            └── {sample}/
                ├── {sample}.input.{chr}.anno.{ref_version}.txt
                └── {sample}.anno.{chr}.{ref_version}_multianno.txt
```

## Notes

- Supports both single-sample and multi-sample analysis
- Can generate GVCF output for joint calling
- Includes PCR indel model handling
- Performs basic annotation using ANNOVAR 