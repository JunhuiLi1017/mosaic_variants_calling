# DeepVariant Workflow Documentation

This document describes the DeepVariant workflow implemented in `rules/variant_calling_deepvariant.smk`. The workflow uses Google's DeepVariant for variant calling.

## Overview

The workflow processes BAM files to identify variants using DeepVariant, a deep learning-based variant caller that can detect SNPs and indels with high accuracy.

## Rules

### 1. run_deepvariant
- **Purpose**: Performs variant calling using DeepVariant
- **Input**: 
  - Sorted and indexed BAM file
  - BAM index file
- **Output**: 
  - VCF file with variants
  - gVCF file with genotype likelihoods
- **Key Parameters**:
  - Model type: WGS
  - Reference genome
  - Chromosome-specific intervals
  - Intermediate results directory
- **Container**: Uses DeepVariant 1.8.0 container

## Dependencies

The workflow requires:
- DeepVariant 1.8.0 (containerized)
- Sorted and indexed BAM files
- Reference genome

## Configuration

The workflow uses several configuration parameters:
- Reference genome
- Interval lists for chromosome-specific processing
- Resource allocation settings
- Model type (WGS)

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
└── 03_variants/
    └── 04_deepvariant/
        └── {sample}/
            ├── {sample}.{chr}.{ref_version}.output.vcf.gz
            └── {sample}.{chr}.{ref_version}.output.gvcf.gz
```

## Notes

- The BAM file must be sorted and indexed
- Duplicate marking is optional (minimal impact on accuracy except at low coverage)
- BQSR is not recommended for DeepVariant 