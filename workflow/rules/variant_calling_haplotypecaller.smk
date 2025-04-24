rule haplotypecaller:
	input:
		bam="{outpath}/02_map/bqsr/{sample}/{sample}.bam"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.haplotypecaller.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		dpsnp138=config['dpsnp138'],
		interval=intervals_dir + "/{chr}.intervals.list",
		pcr_indel_model_command="-pcr_indel_model NONE" if config["pcr_based"] == True else "",
		erc=lambda wildcards: "--emit-ref-confidence NONE" if len(sample_num) == 1 else "--emit-ref-confidence GVCF",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
		HaplotypeCaller \
		-O {output.vcf} \
		-R {params.ref} \
		-I {input.bam}\
		-L {params.interval} \
		{params.pcr_indel_model_command} \
		{params.erc} \
		--dbsnp {params.dpsnp138} > {log} 2>&1
		conda deactivate
		"""

rule mergevcfs:
	input:
		vcf=lambda wildcards: expand(
			"{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz", 
			outpath=wildcards.outpath, 
			sample=wildcards.sample, 
			chr=chromosomes
		),
		tbi=lambda wildcards: expand(
			"{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz.tbi", 
			outpath=wildcards.outpath, 
			sample=wildcards.sample, 
			chr=chromosomes
		)
	output:
		vcf_gz="{outpath}/03_variants/02_haplotypecaller/{sample}/{sample}.vcf.gz",
		vcf_tbi="{outpath}/03_variants/02_haplotypecaller/{sample}/{sample}.vcf.gz.tbi"
	params:
		gatk=config['gatk_current_using'],
		input_args=lambda wildcards, input: " ".join(f"--INPUT {vcf}" for vcf in input.vcf),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
		MergeVcfs \
		{params.input_args} \
		--OUTPUT {output.vcf_gz}
		conda deactivate
		'''

rule haplotypecaller_anno:
	input:
		vcf_gz="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz",
		vcf_tbi="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz.tbi"
	output:
		input_anno="{outpath}/03_variants/02_haplotypecaller/sub_anno/{sample}/{sample}.input.{chr}.anno.{ref_version}.txt",
		vcf_anno="{outpath}/03_variants/02_haplotypecaller/sub_anno/{sample}/{sample}.anno.{chr}.{ref_version}_multianno.txt"
	params:
		ref_version=config["ref_version"],
		data_anno="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar",
		annovar="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl",
		outdir="{outpath}/03_variants/02_haplotypecaller/sub_anno/{sample}/",
		outputanno="{outpath}/03_variants/02_haplotypecaller/sub_anno/{sample}/{sample}.anno.{chr}"
	log:
		"{outpath}/logs/haplotypecaller_anno/{sample}.{chr}.{ref_version}.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		'''
		zcat {input.vcf_gz} | grep -v "#" | awk '{{ print $1"\\t"$2"\\t"$2+length($4)-1"\\t"$4"\\t"$5"\\t"$3 }}' > {output.input_anno}
		perl {params.annovar} \
		{output.input_anno} \
		{params.data_anno}/humandb_{ref_version} \
		-buildver {ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol refGene,dbnsfp42a \
		-operation g,f \
		-nastring . 
		'''

