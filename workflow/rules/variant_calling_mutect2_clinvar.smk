rule mutect2_filter:
	input:
		"{outpath}/03_variants/01_mutect2/result/{sample}/{sample}.mt2pon.filter.vcf.gz"
	output:
		vcf="{outpath}/03_variants/01_mutect2/result/00_somatic_anno/01_filter/{sample}.{ref_version}.pass.output.vcf.gz"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.mutect2_filter.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools view -f PASS {input} -o {output.vcf} -Oz
		tabix -p vcf {output.vcf}
		conda deactivate
		'''

rule annotate_clinvar:
	input:
		vcf="{outpath}/03_variants/01_mutect2/result/00_somatic_anno/01_filter/{sample}.{ref_version}.pass.output.vcf.gz"
	output:
		txt="{outpath}/03_variants/01_mutect2/result/00_somatic_anno/02_annovar_clinvar/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.m2.annotate_clinvar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/01_mutect2/result/00_somatic_anno/02_annovar_clinvar/{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		"""
		perl {params.annovar_dir}/table_annovar.pl \
		{input.vcf} \
		{params.annovar_dir}/humandb_{ref_version} \
		-buildver {ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol refGene,dbnsfp42a,clinvar_20240917,gnomad30_genome,collins_dosage \
		-operation g,f,f,f,f \
		-nastring . \
		-vcfinput
		"""