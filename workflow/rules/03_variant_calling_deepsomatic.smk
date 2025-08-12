rule run_deepsomatic:
	message:
		"The BAM file must be also sorted and indexed. Duplicate marking may be performed, in our analyses there is almost no difference in accuracy except at lower (<20x) coverages. Finally, we recommend that you do not perform BQSR."
	input:
		bam="{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam",
		bai="{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam.bai"
	output:
		vcf="{outpath}/03_variants/deepsomatic/01_raw/{sample}.{chr}.vcf.gz",
		tbi="{outpath}/03_variants/deepsomatic/01_raw/{sample}.{chr}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.deepsomatic.log"
	params:
		model_type="WGS_TUMOR_ONLY",
		ref=config['reference'],
		interval_list=intervals_dir + "/{chr}.intervals.list",
		intermediate_dir="{outpath}/03_variants/deepsomatic/01_raw/{sample}.{chr}.intermediate_results_dir",
		sample="{sample}",
		chr="{chr}"
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	container:
		config["deepsomatic_1.8.0"]
	shell:
		"""
		run_deepsomatic \
		  --model_type={params.model_type} \
		  --ref={params.ref} \
		  --reads_tumor={input.bam} \
		  --output_vcf={output.vcf} \
		  --sample_name_tumor={params.sample} \
		  --num_shards={threads} \
		  --logging_dir={log} \
		  --use_default_pon_filtering=true \
		  --intermediate_results_dir={params.intermediate_dir} \
		  --regions={params.chr}
		"""
