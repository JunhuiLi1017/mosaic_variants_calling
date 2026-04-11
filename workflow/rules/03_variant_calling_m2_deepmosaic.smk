rule pass_filter:
	input:
		vcf="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz",
		index="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/deepmosaic/01_pass_filter/{sample}.mt2pon.pass.vcf.gz"
	log:
		log="{outpath}/03_variants/logs/{sample}.deepmosaic.pass_filter.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image["bcftools_1.9"]
	shell:
		"""
		bcftools view -f PASS {input.vcf} -Oz -o {output.vcf} 2> {log}
		bcftools index {output.vcf}
		"""

rule split_vcf_by_chr:
	input:
		vcf="{outpath}/03_variants/deepmosaic/01_pass_filter/{sample}.mt2pon.pass.vcf.gz"
	output:
		vcf="{outpath}/03_variants/deepmosaic/02_split_chr/{sample}.{chr}.mt2pon.pass.vcf.gz"
	log:
		log="{outpath}/03_variants/logs/{sample}.{chr}.deepmosaic.split_chr.log"
	params:
		interval=intervals_dir + "/{chr}.intervals.list"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image["bcftools_1.9"]
	shell:
		"""
		bcftools view -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf} 2> {log}
		bcftools index {output.vcf}
		"""

rule generate_input_deepmosaic:
	input:
		vcf="{outpath}/03_variants/deepmosaic/02_split_chr/{sample}.{chr}.mt2pon.pass.vcf.gz",
		bam="{outpath}/02_map/08_bqsr/{sample}.bam",
		bai="{outpath}/02_map/08_bqsr/{sample}.bai"
	output:
		txt="{outpath}/03_variants/deepmosaic/03_deepmosaic_input_meta/{sample}.{chr}.deepmosaic_input.txt"
	log:
		log="{outpath}/03_variants/logs/{sample}.{chr}.03_deepmosaic_input_meta.log"
	params:
		depth=lambda wildcards: units.loc[units['sample'] == wildcards.sample, 'depth'].dropna().iloc[0] if len(units.loc[units['sample'] == wildcards.sample, 'depth'].dropna()) > 0 else '200',
		sex=lambda wildcards: units.loc[units['sample'] == wildcards.sample, 'sex'].dropna().iloc[0] if len(units.loc[units['sample'] == wildcards.sample, 'sex'].dropna()) > 0 else 'M'
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		# Write header
		echo -e "#sample_name\\tbam\\tvcf\\tdepth\\tsex" > {output.txt}
		sample_name=$(basename {input.vcf} .mt2pon.pass.vcf.gz)
		echo -e "$sample_name\\t{input.bam}\\t{input.vcf}\\t{params.depth}\\t{params.sex}" >> {output.txt}
		"""

rule deepmosaic_draw:
	input:
		txt="{outpath}/03_variants/deepmosaic/03_deepmosaic_input_meta/{sample}.{chr}.deepmosaic_input.txt"
	output:
		txt="{outpath}/03_variants/deepmosaic/04_deepmosaic_draw/{sample}.{chr}/features.txt"
	log:
		log="{outpath}/03_variants/logs/{sample}.{chr}.04_deepmosaic_draw.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		output_dir="{outpath}/03_variants/deepmosaic/04_deepmosaic_draw/{sample}.{chr}"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image['deepmosaic_1.0.0']
	shell:
		"""
		deepmosaic-draw -i {input.txt} -o {params.output_dir} -a {params.annovar_dir} -b {params.ref_version} -db gnomad41_genome 2> {log}
		"""

rule deepmosaic_predict:
	input:
		txt="{outpath}/03_variants/deepmosaic/04_deepmosaic_draw/{sample}.{chr}/features.txt"
	output:
		"{outpath}/03_variants/deepmosaic/05_deepmosaic_predict/{sample}.{chr}.output.txt"
	log:
		log="{outpath}/03_variants/logs/{sample}.{chr}.05_deepmosaic_predict.log"
	params:
		ref_version=config['ref_version']
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image['deepmosaic_1.0.0']
	shell:
		"""
		deepmosaic-predict -i {input.txt} -o {output} -gb {params.ref_version} 2> {log}
		"""

rule deepmosaic_merge:
	input:
		txt=expand("{{outpath}}/03_variants/deepmosaic/05_deepmosaic_predict/{{sample}}.{chr}.output.txt", chr=CHROMOSOMES)
	output:
		"{outpath}/03_variants/deepmosaic/06_deepmosaic_merge/{sample}.output.txt"
	log:
		log="{outpath}/03_variants/logs/{sample}.06_deepmosaic_merge.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		cat <(cat {input.txt} | head -n 1) <(awk '{{if($22=="mosaic"){{print $0}}}}' {input.txt}) > {output}
		"""

rule deepmosaic_subvcf:
	input:
		txt="{outpath}/03_variants/deepmosaic/06_deepmosaic_merge/{sample}.output.txt",
		vcf="{outpath}/03_variants/deepmosaic/01_pass_filter/{sample}.mt2pon.pass.vcf.gz"
	output:
		"{outpath}/03_variants/deepmosaic/07_deepmosaic_vcf/{sample}.deepmosaic.vcf.gz"
	log:
		log="{outpath}/03_variants/logs/{sample}.07_deepmosaic_vcf.log"
	params:
		bed="{outpath}/03_variants/deepmosaic/07_deepmosaic_vcf/{sample}.deepmosaic.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image["bcftools_1.9"]
	shell:
		"""
		grep -v "#" {input.txt} | awk '{{print $3"\\t"$4-1"\\t"$4}}' > {params.bed}
		bcftools view -R {params.bed} -v snps -Oz -o {output} {input.vcf}
		"""

rule deepmosaic_anno:
	input:
		vcf="{outpath}/03_variants/deepmosaic/07_deepmosaic_vcf/{sample}.deepmosaic.vcf.gz"
	output:
		"{outpath}/03_variants/deepmosaic/08_deepmosaic_anno/{sample}.{ref_version}_multianno.txt"
	log:
		log="{outpath}/03_variants/logs/{sample}.{ref_version}.08_deepmosaic_anno.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/deepmosaic/08_deepmosaic_anno/{sample}"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image["terra_perl_anno"]
	shell:
		"""
		perl {params.annovar_dir}/table_annovar.pl \
		{input.vcf} \
		{params.annovar_dir}/humandb_{params.ref_version} \
		-buildver {params.ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol refGene,dbnsfp42a,clinvar_20240917,gnomad41_genome,gnomad41_exome \
		-operation g,f,f,f,f \
		-nastring . \
		-vcfinput
		"""