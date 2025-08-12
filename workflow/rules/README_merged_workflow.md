# Merged MosaicForecast Workflow

This merged Snakemake workflow combines the Refine and Phase model workflows into a single file, allowing you to run either model based on the `model` parameter.

## Overview

The merged workflow (`03_variant_calling_m2_mosaicforecast_anno_merged.smk`) supports both:
- **Refine model**: Uses `refine_beta` configuration and filters for `mosaic` variants
- **Phase model**: Uses `phase_model` configuration and filters for `hap=3` variants

## Key Differences Between Models

### Refine Model
- Uses `config['refine_beta']` for prediction
- Filters SNV predictions for `$35=="mosaic"`
- Filters INS/DEL predictions for `$35=="hap=3"`
- Output paths: `{outpath}/03_variants/01_mutect2/02_mosaicforecast/...`

### Phase Model
- Uses `config['phase_model']` for prediction
- Filters all predictions (SNV/INS/DEL) for `$35=="hap=3"`
- Output paths: `{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/...`

## Usage

### Running with Refine Model
```bash
snakemake --snakefile workflow/rules/03_variant_calling_m2_mosaicforecast_anno_merged.smk \
  --config model=Refine \
  --configfile config.yaml \
  --cores 8
```

### Running with Phase Model
```bash
snakemake --snakefile workflow/rules/03_variant_calling_m2_mosaicforecast_anno_merged.smk \
  --config model=Phase \
  --configfile config.yaml \
  --cores 8
```

## Required Configuration

Make sure your `config.yaml` includes these parameters:

```yaml
# For Refine model
refine_beta: "/path/to/refine_beta.Rds"
prediction_indel_model: "/path/to/indel_model.Rds"

# For Phase model
phase_model: "/path/to/phase_model.Rds"

# Common parameters
prediction_r_script: "/path/to/Prediction.R"
process_tier_script: "/path/to/process_tier.py"
reformat_script: "/path/to/reformat.py"
```

## Output Structure

The workflow creates model-specific output directories:

### Refine Model Outputs
```
{outpath}/03_variants/01_mutect2/02_mosaicforecast/
├── 07_mosaic_prediction/Refine/
├── 08_anno_gnomAD_dbsnp/Refine/
├── 09_mosaic_tier/Refine/
├── 10_annovar/Refine/
├── 11_annovar_reformat/Refine/
└── Refine/
    ├── {sample}.{cov}.SNV.{ref_version}.mosaic_summary.txt
    ├── {sample}.{cov}.INS.{ref_version}.mosaic_summary.txt
    └── {sample}.{cov}.DEL.{ref_version}.mosaic_summary.txt
```

### Phase Model Outputs
```
{outpath}/03_variants/01_mutect2/02_mosaicforecast/
├── 07_mosaic_prediction/Phase/
├── 08_anno_gnomAD_dbsnp/Phase/
├── 09_mosaic_tier/Phase/
├── 10_annovar/Phase/
├── 11_annovar_reformat/Phase/
└── Phase/
    ├── {sample}.{cov}.SNV.{ref_version}.mosaic_summary.txt
    ├── {sample}.{cov}.INS.{ref_version}.mosaic_summary.txt
    └── {sample}.{cov}.DEL.{ref_version}.mosaic_summary.txt
```

## Key Features

1. **Conditional Logic**: The workflow uses shell conditionals to choose between Refine and Phase model parameters and filtering criteria.

2. **Model-Specific Paths**: Output directories include the model name to keep results separate.

3. **Unified Interface**: Single workflow file handles both models, reducing code duplication.

4. **Backward Compatibility**: Maintains the same output structure as the original separate workflows.

## Migration from Separate Workflows

If you were previously using separate workflow files:

1. **Old Refine workflow**: `03_variant_calling_m2_mosaicforecast_anno.smk`
2. **Old Phase workflow**: `variant_calling_m2_mf_phase.smk`

Replace with:
```bash
# Instead of running separate workflows, use:
snakemake --snakefile workflow/rules/03_variant_calling_m2_mosaicforecast_anno_merged.smk --config model=Refine
snakemake --snakefile workflow/rules/03_variant_calling_m2_mosaicforecast_anno_merged.smk --config model=Phase
```

## Troubleshooting

### Common Issues

1. **Model parameter not set**: Ensure you specify `model=Refine` or `model=Phase` in your command or config.

2. **Missing configuration**: Verify all required model-specific parameters are defined in your config file.

3. **Output directory conflicts**: The workflow creates separate directories for each model, so there should be no conflicts.

### Validation

To verify the workflow is working correctly:

1. Check that the correct model-specific parameters are being used in the logs
2. Verify that output files are created in the expected model-specific directories
3. Confirm that the filtering criteria match the expected model behavior

## Dependencies

The merged workflow requires the same dependencies as the original workflows:
- GATK4
- bcftools
- bedtools
- ANNOVAR
- R (for prediction scripts)
- Python (for processing scripts) 