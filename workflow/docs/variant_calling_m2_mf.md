# Mutect2 Workflow Documentation

This document describes the Mutect2 workflow implemented in `rules/variant_calling_m2_mf.smk`. The workflow uses GATK's Mutect2 for somatic variant calling.

## Overview

The workflow processes BAM files to identify somatic variants using Mutect2, followed by filtering and annotation steps. It is specifically designed for detecting mosaic variants.

## Rules

### 1. Mutect2
- **Purpose**: Performs somatic variant calling using Mutect2
- **Input**: 
  - BAM file
  - Panel of normals (PON)
  - Germline resource
- **Output**: 
  - VCF file with raw variants
  - Statistics file
- **Key Parameters**:
  - Reference genome
  - PON file
  - Germline resource
  - Chromosome-specific intervals
- **Tools**: GATK 4.6.1.0

### 2. vcf_merge
- **Purpose**: Merges chromosome-specific VCF files
- **Input**: Individual chromosome VCF files
- **Output**: Merged VCF file
- **Tools**: BCFtools

### 3. merge_mutectstats
- **Purpose**: Merges statistics files from individual chromosomes
- **Input**: Individual chromosome statistics files
- **Output**: Merged statistics file
- **Tools**: Custom Python script

### 4. FilterMutectCall
- **Purpose**: Applies filtering to the merged VCF
- **Input**: 
  - Merged VCF file
  - Statistics file
- **Output**: Filtered VCF file
- **Tools**: GATK FilterMutectCalls

### 5. MT2_initial_filter
- **Purpose**: Initial filtering of Mutect2 calls
- **Input**: Filtered VCF file
- **Output**: BED file with filtered variants
- **Tools**: Custom Python script

### 6. repeat_filter
- **Purpose**: Filters out variants in segmental duplications
- **Input**: BED file
- **Output**: Filtered BED file
- **Tools**: BEDtools

### 7. annovar_formatter
- **Purpose**: Formats variants for ANNOVAR annotation
- **Input**: BED file
- **Output**: Formatted file for ANNOVAR
- **Tools**: Custom Python script

### 8. MAF0_extraction Rules
- **Purpose**: Extracts different variant types
- **Input**: Formatted file
- **Output**: BED files for SNVs, insertions, and deletions
- **Tools**: Custom Python script

### 9. feature_extraction Rules
- **Purpose**: Extracts features for variant classification
- **Input**: 
  - BED files
  - BAM files
- **Output**: Feature files
- **Tools**: Custom Python scripts

### 10. g1000_avail_acess_filter Rules
- **Purpose**: Filters variants based on 1000 Genomes accessibility
- **Input**: Feature files
- **Output**: Filtered feature files
- **Tools**: Custom Python script

### 11. Prediction Rules
- **Purpose**: Predicts mosaic variants
- **Input**: Filtered feature files
- **Output**: Prediction files
- **Tools**: Custom Python script

### 12. extract_bed Rules
- **Purpose**: Extracts predicted variants to BED format
- **Input**: Prediction files
- **Output**: BED files
- **Tools**: Custom Python script

### 13. extrac_subvcf Rules
- **Purpose**: Extracts variants to VCF format
- **Input**: 
  - BED files
  - Original VCF file
- **Output**: VCF files
- **Tools**: BCFtools

### 14. anno_gnomAD_dbsbp Rules
- **Purpose**: Annotates variants with gnomAD and dbSNP information
- **Input**: VCF files
- **Output**: Annotated VCF files
- **Tools**: BCFtools

### 15. process_tier
- **Purpose**: Processes variant tiers
- **Input**: Annotated VCF files
- **Output**: Tiered VCF files
- **Tools**: Custom Python script

### 16. Annotation Rules
- **Purpose**: Performs ANNOVAR annotation
- **Input**: Tiered VCF files
- **Output**: Annotated files
- **Tools**: ANNOVAR

### 17. reformat_annotation_clinvar
- **Purpose**: Reformats annotation with ClinVar information
- **Input**: Annotated files
- **Output**: Reformatted files
- **Tools**: Custom Python script

### 18. annotate_rcnv_gnomadlof_snv
- **Purpose**: Annotates SNVs with rCNV and gnomAD LoF information
- **Input**: Reformatted files
- **Output**: Annotated files
- **Tools**: Custom Python script

### 19. reformat_rcnv_gnomadlof
- **Purpose**: Reformats rCNV and gnomAD LoF annotations
- **Input**: Annotated files
- **Output**: Reformatted files
- **Tools**: Custom Python script

### 20. summary
- **Purpose**: Generates summary statistics
- **Input**: Various annotation files
- **Output**: Summary files
- **Tools**: Custom Python script

## Dependencies

The workflow requires several external tools:
- GATK 4.6.1.0
- ANNOVAR
- BEDtools
- BCFtools
- Python 3.7.1
- R

## Configuration

The workflow uses several configuration parameters:
- Reference genome
- Panel of normals
- Germline resource
- Resource allocation settings
- Model parameters for variant classification

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
└── 03_variants/
    └── 01_mutect2/
        ├── 01_somatic/
        │   ├── 01_result/
        │   │   ├── sub/
        │   │   │   └── {sample}.{chr}.mt2pon.vcf.gz
        │   │   └── {sample}.mt2pon.merged.vcf.gz
        │   ├── 02_filtered/
        │   ├── 03_annovar/
        │   ├── 04_feature/
        │   ├── 05_prediction/
        │   ├── 06_bed/
        │   ├── 07_vcf/
        │   ├── 08_anno/
        │   ├── 09_tier/
        │   ├── 10_annovar/
        │   ├── 11_reformat/
        │   ├── 12_rcnv_gnomadlof/
        │   └── 13_summary/
        └── 02_germline/
            └── ...
```

## Notes

- The workflow is specifically designed for detecting mosaic variants
- Includes extensive filtering and annotation steps
- Uses machine learning for variant classification
- Generates comprehensive summary statistics 