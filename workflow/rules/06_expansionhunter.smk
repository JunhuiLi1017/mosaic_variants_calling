rule expansionhunter_v5_0_0:
	input:
		bam="{outpath}/02_map/08_bqsr/{sample}.bam",
		bai="{outpath}/02_map/08_bqsr/{sample}.bam.bai"
	output:
		vcf="{outpath}/05_expansionhunter/01_expansionhunter/{sample}.vcf"
	log:
		"{outpath}/05_expansionhunter/logs/{sample}.expansionhunter.log"
	params:
		ref=config['reference'],
		variant_catalog=config['variant_catalog'],
        sex=[units[units['sample'] == sample]['sex'].iloc[0] for sample in Sample.unique()],
		output_prefix="{outpath}/05_expansionhunter/01_expansionhunter/{sample}"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image["expansionhunter_5.0.0"]
	shell:
		"""
		ExpansionHunter --reads {input.bam} \
				--reference {params.ref} \
                --sex {params.sex} \
				--variant-catalog {params.variant_catalog} \
				--output-prefix {params.output_prefix}
		"""

rule stranger:
	input:
		vcf="{outpath}/05_expansionhunter/01_expansionhunter/{sample}.vcf",
	output:
		vcf="{outpath}/05_expansionhunter/02_stranger/{sample}.stranger.vcf"
	log:
		"{outpath}/05_expansionhunter/logs/{sample}.stranger.log"
	params:
		ref=config['reference'],
		variant_catalog=config['variant_catalog']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image["stranger_0.10.2"]
	shell:
		"""
		stranger --vcf {input.vcf} \
				--output {output.vcf} \
				--repeats-file {params.variant_catalog}
		"""