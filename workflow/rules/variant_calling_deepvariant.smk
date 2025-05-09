rule run_deepvariant:
	message:
		"The BAM file must be also sorted and indexed. Duplicate marking may be performed, in our analyses there is almost no difference in accuracy except at lower (<20x) coverages. Finally, we recommend that you do not perform BQSR."
	input:
		bam="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam",
		bai="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam.bai"
	output:
		vcf="{outpath}/03_variants/04_deepvariant/{sample}/{sample}.{individual_chr}.{ref_version}.output.vcf.gz",
		gvcf="{outpath}/03_variants/04_deepvariant/{sample}/{sample}.{individual_chr}.{ref_version}.output.gvcf.gz"
	log:
		"{outpath}/03_variants/logs/{sample}.{individual_chr}.{ref_version}.deepvariant.log"
	params:
		model_type="WGS",
		ref=config['reference'],
		interval_list=intervals_dir + "/{individual_chr}.intervals.list",
		intermediate_dir="{outpath}/03_variants/04_deepvariant/{sample}/{sample}.{individual_chr}.{ref_version}.intermediate_results_dir",
		sample="{sample}",
		chr="{individual_chr}"
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	container:
		"library://junhuili/deepvariant/deepvariant:1.8.0"
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
