# DeepSomatic Workflow Documentation

This document describes the DeepSomatic workflow implemented in `rules/variant_calling_deepsomatic.smk`. The workflow uses DeepSomatic for somatic variant calling.

## Overview

The workflow processes BAM files to identify somatic variants using DeepSomatic, a deep learning-based variant caller specifically designed for somatic variant detection.

## Rules

### 1. run_deepsomatic
- **Purpose**: Performs somatic variant calling using DeepSomatic
- **Input**: 
  - Sorted and indexed BAM file
  - BAM index file
- **Output**: 
  - VCF file with somatic variants
- **Key Parameters**:
  - Model type: WGS_TUMOR_ONLY
  - Reference genome
  - Chromosome-specific intervals
  - Intermediate results directory
- **Container**: Uses DeepSomatic 1.8.0 container

## Dependencies

The workflow requires:
- DeepSomatic 1.8.0 (containerized)
- Sorted and indexed BAM files
- Reference genome

## Configuration

The workflow uses several configuration parameters:
- Reference genome
- Interval lists for chromosome-specific processing
- Resource allocation settings
- Model type (WGS_TUMOR_ONLY)

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
└── 03_variants/
    └── 03_deepsomatic/
        └── {sample}/
            └── {sample}.{chr}.{ref_version}.output.vcf.gz
```

## Notes

- The BAM file must be sorted and indexed
- Duplicate marking is optional (minimal impact on accuracy except at low coverage)
- BQSR is not recommended for DeepSomatic
- Uses default PON (Panel of Normals) filtering 