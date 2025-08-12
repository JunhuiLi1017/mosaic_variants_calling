# Mutect2 PON Workflow Documentation

This document describes the Mutect2 PON (Panel of Normals) workflow implemented in `rules/variant_calling_m2_mf_pon.smk`. The workflow creates a panel of normals for Mutect2 somatic variant calling.

## Overview

The workflow processes normal sample BAM files to create a panel of normals (PON) for Mutect2 somatic variant calling. The PON helps filter out common sequencing artifacts and germline variants.

## Rules

### 1. GetPileupSummaries
- **Purpose**: Generates pileup summaries for normal samples
- **Input**: 
  - BAM file
  - BAM index file
- **Output**: 
  - Pileup summary file
  - Index file
- **Key Parameters**:
  - Reference genome
  - Variant sites file
- **Tools**: GATK 4.6.1.0

### 2. CalculateContamination
- **Purpose**: Calculates contamination estimates
- **Input**: Pileup summary file
- **Output**: Contamination table
- **Tools**: GATK CalculateContamination

### 3. Mutect2_pon
- **Purpose**: Creates panel of normals
- **Input**: 
  - BAM files from normal samples
  - Contamination tables
- **Output**: 
  - PON VCF file
  - Index file
- **Key Parameters**:
  - Reference genome
  - Germline resource
  - Chromosome-specific intervals
- **Tools**: GATK 4.6.1.0

## Dependencies

The workflow requires:
- GATK 4.6.1.0
- Reference genome
- Germline resource (gnomAD)
- Normal sample BAM files

## Configuration

The workflow uses several configuration parameters:
- Reference genome
- Germline resource
- Resource allocation settings
- Chromosome-specific intervals

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
└── 03_variants/
    └── 01_mutect2/
        └── 00_pon/
            ├── pileup_summaries/
            │   └── {sample}/
            │       └── {sample}.{chr}.getpileupsummaries.table
            ├── contamination/
            │   └── {sample}/
            │       └── {sample}.{chr}.calculatecontamination.table
            └── pon/
                └── {chr}/
                    └── {chr}.pon.vcf.gz
```

## Notes

- The PON is created from normal samples only
- Contamination estimation is performed for each sample
- The PON is created per chromosome
- The PON can be used in subsequent Mutect2 somatic variant calling 