rule run_deepvariant:
	message:
		"The BAM file must be also sorted and indexed. Duplicate marking may be performed, in our analyses there is almost no difference in accuracy except at lower (<20x) coverages. Finally, we recommend that you do not perform BQSR."
	input:
		bam="{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam",
		bai="{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam.bai"
	output:
		vcf="{outpath}/03_variants/deepvariant/01_raw/{sample}.{chr}.vcf.gz",
		gvcf="{outpath}/03_variants/deepvariant/01_raw/{sample}.{chr}.gvcf.gz",
		tbi="{outpath}/03_variants/deepvariant/01_raw/{sample}.{chr}.vcf.gz.tbi",
		gtbi="{outpath}/03_variants/deepvariant/01_raw/{sample}.{chr}.gvcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.deepvariant.log"
	params:
		model_type="WGS",
		ref=config['reference'],
		interval_list=intervals_dir + "/{chr}.intervals.list",
		intermediate_dir="{outpath}/03_variants/deepvariant/01_raw_sub_vcf/{sample}.{chr}.intermediate_results_dir",
		sample="{sample}",
		chr="{chr}"
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	container:
		"docker://ksinghal28/deepvariant:1.8.1"
	shell:
		"""
		run_deepvariant \
		  --model_type={params.model_type} \
		  --ref={params.ref} \
		  --reads={input.bam} \
		  --output_vcf={output.vcf} \
		  --output_gvcf={output.gvcf} \
		  --num_shards={threads} \
		  --logging_dir={log} \
		  --intermediate_results_dir={params.intermediate_dir} \
		  --regions={params.chr}
		"""