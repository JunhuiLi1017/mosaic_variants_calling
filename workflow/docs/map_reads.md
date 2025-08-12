# Read Mapping Workflow Documentation

This document describes the read mapping workflow implemented in `rules/map_reads.smk`. The workflow handles the alignment of sequencing reads to a reference genome using BWA-MEM and includes subsequent processing steps.

## Overview

The workflow processes FASTQ files to align reads to a reference genome, followed by sorting, duplicate marking, and base quality score recalibration (BQSR).

## Rules

### 1. map_reads
- **Purpose**: Aligns reads to reference genome using BWA-MEM
- **Input**: FASTQ files (R1 and R2 for paired-end)
- **Output**: 
  - Raw BAM file
- **Key Parameters**:
  - Read group information
  - Reference genome
  - Thread count based on resource allocation

### 2. map_reads_sort
- **Purpose**: Sorts aligned reads by coordinate
- **Input**: Raw BAM file
- **Output**: 
  - Sorted BAM file
  - Statistics file
- **Tools**: Sambamba

### 3. map_reads_frag
- **Purpose**: Analyzes fragment size distribution
- **Input**: Sorted BAM file
- **Output**: 
  - Fragment size plot (PNG)
  - Fragment size statistics (TXT)
- **Tools**: DeepTools

### 4. remove_dup
- **Purpose**: Marks or removes duplicate reads
- **Input**: Sorted BAM files
- **Output**: 
  - Deduplicated BAM file
  - Index file
  - Duplicate metrics
- **Parameters**:
  - PCR-based or non-PCR-based duplicate handling
  - Optical duplicate distance

### 5. SetNmMdAndUqTags
- **Purpose**: Sets NM, MD, and UQ tags in BAM file
- **Input**: Deduplicated BAM file
- **Output**: Tagged BAM file
- **Tools**: GATK

### 6. sortbam
- **Purpose**: Final sorting of BAM file
- **Input**: Tagged BAM file
- **Output**: 
  - Sorted BAM file
  - Index file
- **Tools**: Sambamba

### 7. BaseRecalibrator
- **Purpose**: Generates base quality score recalibration tables
- **Input**: Sorted BAM file
- **Output**: Recalibration tables per chromosome
- **Parameters**:
  - Known SNPs (dbSNP)
  - Known indels
  - Mills and 1000G indels

### 8. GatherBQSRreports
- **Purpose**: Combines recalibration tables
- **Input**: Individual chromosome recalibration tables
- **Output**: Combined recalibration table

### 9. ApplyBQSR
- **Purpose**: Applies base quality score recalibration
- **Input**: 
  - Sorted BAM file
  - Recalibration table
- **Output**: 
  - Recalibrated BAM file per chromosome
  - Statistics file

### 10. GatherBQSRBam
- **Purpose**: Combines recalibrated BAM files
- **Input**: Individual chromosome BAM files
- **Output**: 
  - Combined BAM file
  - Index file

## Dependencies

The workflow requires several external tools:
- BWA-MEM
- Sambamba
- GATK 4.6.1.0
- DeepTools
- SAMtools

## Configuration

The workflow uses several configuration parameters:
- Reference genome
- Known variant databases (dbSNP, Mills, 1000G)
- Resource allocation settings
- PCR-based vs non-PCR-based settings

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
├── 02_map/
│   ├── bwa/
│   │   ├── raw/
│   │   └── sort/
│   ├── dup/
│   ├── settags/
│   └── bqsr/
│       ├── sub_recal/
│       ├── sub_bqsr/
│       └── {sample}/
└── logs/
    └── bwa/
``` 