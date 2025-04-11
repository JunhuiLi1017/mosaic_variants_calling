rule all:
	input:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/02ProcessedData/WGS/01_Mosaic_variants/test_deep_somatic/output1/sompy_output1/output1.vcf.gz"

rule run_deepsomatic:
	message:
		"The BAM file must be also sorted and indexed. Duplicate marking may be performed, in our analyses there is almost no difference in accuracy except at lower (<20x) coverages. Finally, we recommend that you do not perform BQSR."
	input:
		bam="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/02ProcessedData/WGS/01_Mosaic_variants/test_deep_somatic/input/data/HCC1395_wgs.tumor.chr1.bam",
		bai="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/02ProcessedData/WGS/01_Mosaic_variants/test_deep_somatic/input/data/HCC1395_wgs.tumor.chr1.bam.bai"
	output:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/02ProcessedData/WGS/01_Mosaic_variants/test_deep_somatic/output1/sompy_output1/output1.vcf.gz"
	log:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/02ProcessedData/WGS/01_Mosaic_variants/test_deep_somatic/output1/sompy_output1/log"
	params:
		model_type="WGS_TUMOR_ONLY",
		ref="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/02ProcessedData/WGS/01_Mosaic_variants/test_deep_somatic/input/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna",
		intermediate_dir="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/02ProcessedData/WGS/01_Mosaic_variants/test_deep_somatic/output1/sompy_output1/tmp_dir"
	threads:
		8
	resources:
		mem_mb=2000
	container:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/mosaic_variants_calling_largeinput/workflow/envs/deepsomatic_1.8.0.sif"
	shell:
		"""
		run_deepsomatic \
		  --model_type={params.model_type} \
		  --ref={params.ref} \
		  --reads_tumor={input.bam} \
		  --output_vcf={output} \
		  --sample_name_tumor="HCC1395Tumor" \
		  --num_shards={threads} \
		  --logging_dir={log} \
		  --use_default_pon_filtering=true \
		  --intermediate_results_dir={params.intermediate_dir} \
		  --regions=chr1
		"""