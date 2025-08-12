# Quality Control Workflow Documentation

This document describes the quality control workflow implemented in `rules/qc.smk`. The workflow performs quality control checks on sequencing data and generates comprehensive reports.

## Overview

The workflow processes raw sequencing data through multiple QC tools to assess data quality and generate reports. It includes read trimming, quality assessment, and multi-sample report generation.

## Rules

### 1. fastp
- **Purpose**: Performs read trimming and quality control
- **Input**: Raw FASTQ files
- **Output**: 
  - Trimmed FASTQ files
  - JSON report
  - HTML report
- **Parameters**:
  - Trim settings for R1 and R2
  - Thread count based on resource allocation

### 2. fastqc
- **Purpose**: Performs quality control on sequencing data
- **Input**: Trimmed FASTQ files
- **Output**: 
  - HTML reports
  - ZIP archives with detailed QC metrics
- **Tools**: FastQC

### 3. samtools_stats
- **Purpose**: Generates statistics for BAM files
- **Input**: BAM file after BQSR
- **Output**: Statistics file
- **Tools**: SAMtools

### 4. multiqc
- **Purpose**: Aggregates results from multiple QC tools
- **Input**: Results from various QC tools
- **Output**: 
  - MultiQC HTML report
- **Tools**: MultiQC

## Dependencies

The workflow requires several external tools:
- fastp
- FastQC
- SAMtools
- MultiQC

## Configuration

The workflow uses several configuration parameters:
- Trim settings for R1 and R2 reads
- Resource allocation settings
- Output directory structure

## Output Structure

The workflow generates outputs in the following directory structure:
```
{outpath}/
├── 01_multiqc/
│   ├── fastp/
│   │   ├── {sample}_{library}_{flowlane}.R1.fastq.gz
│   │   ├── {sample}_{library}_{flowlane}.R2.fastq.gz
│   │   ├── {sample}_{library}_{flowlane}.json
│   │   └── {sample}_{library}_{flowlane}.html
│   ├── fastqc/
│   │   ├── {sample}_{library}_{flowlane}.R1_fastqc.html
│   │   ├── {sample}_{library}_{flowlane}.R2_fastqc.html
│   │   ├── {sample}_{library}_{flowlane}.R1_fastqc.zip
│   │   └── {sample}_{library}_{flowlane}.R2_fastqc.zip
│   └── multiqc_report.html
├── 02_map/
│   └── bqsr_stat/
│       └── {sample}.sort.rmdup.bqsr.stat
└── logs/
    ├── fastp/
    ├── fastqc/
    └── multiqc/
``` 