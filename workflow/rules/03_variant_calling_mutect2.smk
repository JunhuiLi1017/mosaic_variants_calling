rule Mutect2:
	input:
		"{outpath}/02_map/08_bqsr/{sample}.bam"
	output:
		raw_vcf="{outpath}/03_variants/mutect2/00_initial_call/01_chr_sub/{sample}.{chr}.vcf.gz",
		raw_stat="{outpath}/03_variants/mutect2/00_initial_call/01_chr_sub/{sample}.{chr}.vcf.gz.stats"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.mutect.log"
	params:
		pon=config['pon'],
		af_only_gnomad=config['af_only_gnomad'],
		ref=config['reference'],
		sample="{sample}",
		interval_list=intervals_dir + "/{chr}.intervals.list",
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		'''
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		Mutect2 \
		-R {params.ref} \
		-I {input} \
		--pon {params.pon} \
		-tumor {params.sample} \
		--germline-resource {params.af_only_gnomad} \
		-L {params.interval_list} \
		--interval-padding 100 \
		-O {output.raw_vcf} > {log} 2>&1
		'''

rule vcf_merge:
	message:
		"""
		---
		merge the vcf file for individual chromosome into single one vcf for each sample
		---
		"""
	input:
		lambda wildcards: [f"{wildcards.outpath}/03_variants/mutect2/00_initial_call/01_chr_sub/{wildcards.sample}.{chr}.vcf.gz" for chr in CHROMOSOMES]
	output:
		args="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.args",
		vcf="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz",
		idx="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.merge.logs"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["vcftools_0.1.16"]
	shell:
		'''
		ls {input} > {output.args}
		vcf-concat -f {output.args}  |  bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		'''

rule merge_mutectstats:
	input:
		lambda wildcards: [f"{wildcards.outpath}/03_variants/mutect2/00_initial_call/01_chr_sub/{wildcards.sample}.{chr}.vcf.gz.stats" for chr in CHROMOSOMES]
	output:
		"{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz.stats"
	params:
		list_para = lambda wildcards: ' '.join([f"--stats {wildcards.outpath}/03_variants/mutect2/00_initial_call/01_chr_sub/{wildcards.sample}.{chr}.mt2pon.vcf.gz.stats" for chr in CHROMOSOMES]),
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	log:
		"{outpath}/03_variants/logs/{sample}.MergeMutectStats.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		'''
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		MergeMutectStats \
		{params.list_para} \
		-O {output} > {log} 2>&1
		'''

rule FilterMutectCall:
	input:
		vcf="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz",
		stats="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz.stats",
		index="{outpath}/03_variants/mutect2/00_initial_call/02_merge/{sample}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/mutect2/00_initial_call/03_filter_mutect_call/{sample}.vcf.gz",
		tbi="{outpath}/03_variants/mutect2/00_initial_call/03_filter_mutect_call/{sample}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.filter.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		'''
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		FilterMutectCalls -R {params.ref} -V {input.vcf} -O {output} > {log} 2>&1
		'''

rule spliet_vcf_m2:
	input:
		vcf="{outpath}/03_variants/mutect2/00_initial_call/03_filter_mutect_call/{sample}.vcf.gz",
		tbi="{outpath}/03_variants/mutect2/00_initial_call/03_filter_mutect_call/{sample}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/mutect2/01_raw/{sample}.{chr}.vcf.gz",
		idx="{outpath}/03_variants/mutect2/01_raw/{sample}.{chr}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/m2_split_vcf.{sample}.{chr}.log"
	params:
		sample="{sample}",
		chr="{chr}"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		"""
		 bcftools view \
            -s {params.sample} \
            -r {params.chr} \
            -o {output.vcf} \
            -Oz \
            {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		"""	