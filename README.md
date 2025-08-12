# Mosaic Variants Calling Analysis Pipeline

A comprehensive Snakemake-based pipeline for Mosaic Variants Calling analysis with advanced variant calling capabilities, including mosaic variant detection using MosaicForecast.

## 🚀 Key Features

### 📊 **Resource Management**
- **Flexible Resource Configuration**: Configurable memory and thread allocation through `config/resources.yaml`
- **Tiered Resource Levels**: 
  - Low: 1GB RAM, 2 threads
  - Medium: 2GB RAM, 4 threads  
  - High: 3GB RAM, 8 threads
  - Very High: 4GB RAM, 12 threads
- **Dynamic Resource Assignment**: Each rule can specify appropriate resource requirements

### 🔧 **Unit-Based Configuration**
- **Sample-Level Customization**: Individual sample parameters in `config/units.tsv`
- **Per-Sample Trimming**: Configurable `trim_front1`, `trim_front2`, `trim_tail1`, `trim_tail2` for each sample
- **PCR Status Tracking**: `pcr_based` field for PCR-based vs non-PCR-based samples
- **Library and Flowlane Management**: Organized sample tracking with library and flowlane identifiers

### 📈 **Comprehensive Reporting**
- **MultiQC Integration**: Automated quality control reports for each pipeline stage
- **Mapping Reports**: Detailed alignment statistics and coverage analysis
- **Variant Reports**: Comprehensive variant calling summaries
- **Mosaic Reports**: Specialized reports for mosaic variant analysis
- **HTML Reports**: Interactive web-based reports for easy visualization

### 🐳 **Container Technology**
- **Singularity Support**: Full containerization with `--use-singularity` flag
- **Docker Containers**: Pre-built containers for all major tools:
  - **QC Tools**: FastQC, FastP, MultiQC, BedTools
  - **Mapping**: BWA, Samtools, Sambamba, DeepTools
  - **Variant Calling**: GATK4, DeepVariant, DeepSomatic, HaplotypeCaller
  - **Analysis**: VCFtools, BCFtools, R, Python
- **Version Control**: Specific container versions for reproducibility

### 🔬 **Advanced Variant Calling**
- **Multiple Callers**: Support for multiple variant calling algorithms:
  - **GATK4 Mutect2**: Industry-standard somatic variant caller
  - **DeepVariant**: Google's deep learning-based caller
  - **DeepSomatic**: Specialized somatic variant detection
  - **HaplotypeCaller**: Germline variant calling
- **MosaicForecast Integration**: Advanced mosaic variant detection with:
  - **Refine Model**: For general mosaic variant prediction
  - **Phase Model**: For phase-specific variant analysis
  - **Machine Learning**: Trained models for SNV and indel prediction

### 🧬 **Reference Genome Support**
- **hg38 Configuration**: Complete reference genome setup
- **Database Integration**: Comprehensive annotation databases:
  - dbSNP, gnomAD, ClinVar
  - Mills and 1000G indel sets
  - Panel of Normals (PON)
  - Segmental duplications and repeats
- **Interval Management**: Chromosome-specific interval lists

### 🔄 **Workflow Management**
- **Modular Design**: Organized into logical workflow stages:
  - `01_qc`: Quality control and preprocessing
  - `02_map`: Read mapping and alignment
  - `03_variant_calling`: Multiple variant calling approaches
  - `04_filter_anno`: Variant filtering and annotation
- **Cluster Integration**: LSF cluster support with resource-aware job submission
- **Latency Management**: Configurable latency wait for file system synchronization

### 📋 **Quality Control**
- **FastQC**: Raw read quality assessment
- **FastP**: Adapter trimming and quality filtering
- **MultiQC**: Aggregated QC reports across all samples
- **Coverage Analysis**: Mean depth calculation and coverage statistics

### 🎯 **Mosaic Variants Calling Optimization**
- **Interval-Based Processing**: Chromosome-specific processing for efficiency
- **Coverage Optimization**: Designed for Mosaic Variants Calling data
- **PCR Artifact Handling**: Specialized handling for PCR-based libraries
- **Custom Trimming**: Sample-specific trimming parameters

## 📁 Pipeline Structure

```
workflow/
├── Snakefile                 # Main workflow file
├── rules/                    # Individual workflow modules
│   ├── 00_common.smk        # Common functions and parameters
│   ├── 01_qc.smk            # Quality control
│   ├── 02_map_reads.smk     # Read mapping
│   ├── 02_map_report.smk    # Mapping reports
│   ├── 03_variant_calling_*.smk  # Variant calling modules
│   └── 04_filter_anno_clinvar.smk # Filtering and annotation
├── envs/                     # Conda environment definitions
├── scripts/                  # Custom analysis scripts
└── docs/                     # Documentation

config/
├── config.hg38.yaml         # Main configuration file
├── resources.yaml           # Resource definitions
└── units.tsv               # Sample-specific parameters
```

## 🚀 Quick Start

1. **Configure Resources**:
   ```bash
   # Edit config/resources.yaml for your cluster
   ```

2. **Set Up Samples**:
   ```bash
   # Edit config/units.tsv with your sample information
   ```

3. **Run Pipeline**:
   ```bash
   # Using Singularity containers
   snakemake --use-singularity --configfile config/config.hg38.yaml
   
   # On cluster with LSF
   bash workflow/submit.sh
   ```

## 📊 Output Reports

The pipeline generates comprehensive reports:
- `{outpath}/01_multiqc/multiqc_report.html` - QC summary
- `{outpath}/02_map/multiqc_report_map.html` - Mapping statistics  
- `{outpath}/03_variants/mosaic_report.html` - Mosaic variant analysis
- `{outpath}/03_variants/variant_report.txt` - Variant calling summary

## 🔧 Configuration

### Resource Configuration
```yaml
# config/resources.yaml
resource:
  low:
    mem_mb: 1000
    threads: 2
  medium:
    mem_mb: 2000
    threads: 4
  high:
    mem_mb: 3000
    threads: 8
  very_high:
    mem_mb: 4000
    threads: 12
```

### Sample Configuration
```tsv
# config/units.tsv
sample	library	flowlane	platform	fq1	fq2	trim_front1	trim_front2	trim_tail1	trim_tail2	pcr_based
sample1	A	B	ILLUMINA	/path/to/R1.fastq.gz	/path/to/R2.fastq.gz	9	9	0	0	Yes
```

## 🐳 Container Support

The pipeline uses pre-built containers for all tools:
- **Docker**: `docker://broadinstitute/gatk:4.6.1.0`
- **Singularity**: Automatic conversion from Docker URIs
- **Custom Containers**: Specialized containers for specific tools

## 📈 Performance Features

- **Parallel Processing**: Chromosome-level parallelization
- **Resource Optimization**: Dynamic memory allocation based on input size
- **Temporary File Management**: Efficient handling of intermediate files
- **Cluster Integration**: LSF job submission with resource specifications

## 🔬 Scientific Features

- **Mosaic Variant Detection**: Advanced ML-based mosaic variant calling
- **Multiple Model Support**: Refine and Phase models for different analysis needs
- **Comprehensive Annotation**: Integration with major variant databases
- **Quality Metrics**: Extensive quality control and validation metrics

## 📝 Notes

1. **Latency Wait**: Set to 13500ms for Mutect2 output handling
2. **Dynamic Inputs**: Uses functions/lambda expressions for conditional input resolution
3. **Resource Config**: Centralized resource management for all rules
4. **Unit-Level Settings**: PCR status and trimming parameters per sample

## 🤝 Contributing

This pipeline is designed for Mosaic Variants Calling analysis with a focus on mosaic variant detection. The modular design allows for easy customization and extension.
