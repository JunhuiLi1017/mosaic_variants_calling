# Mutect2 PON Merge Workflow Documentation

This document describes the Mutect2 PON merge workflow implemented in `rules/variant_calling_m2_mf_pon_merge.smk`. The workflow merges chromosome-specific panel of normals (PON) files into a single PON file.

## Overview

The workflow combines individual chromosome PON files into a single merged PON file that can be used for Mutect2 somatic variant calling across the entire genome.

## Rules

### 1. merge_pon
- **Purpose**: Merges chromosome-specific PON files
- **Input**: 
  - Individual chromosome PON VCF files
  - Index files
- **Output**: 
  - Merged PON VCF file
  - Index file
- **Key Parameters**:
  - Reference genome
  - Chromosome list
- **Tools**: GATK MergeVcfs

## Dependencies

The workflow requires:
- GATK 4.6.1.0
- Reference genome
- Chromosome-specific PON files

## Configuration

The workflow uses several configuration parameters:
- Reference genome
- Resource allocation settings
- Chromosome list

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
└── 03_variants/
    └── 01_mutect2/
        └── 00_pon/
            └── pon_merged/
                └── pon.vcf.gz
```

## Notes

- The merged PON file contains variants from all chromosomes
- The merged PON can be used as input for Mutect2 somatic variant calling
- The workflow assumes that individual chromosome PON files are already created
- The merge process preserves all variant information from individual PON files 