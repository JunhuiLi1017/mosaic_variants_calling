rule haplotypecaller:
	input:
		bam="{outpath}/02_map/08_bqsr/{sample}/{sample}.bam"
	output:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/01_sub_vcf/{sample}.{chr}.g.vcf.gz",
		tbi="{outpath}/03_variants/haplotypecaller/00_initial_call/01_sub_vcf/{sample}.{chr}.g.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.haplotypecaller.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		dbsnp138=config['dbsnp138'],
		interval=intervals_dir + "/{chr}.intervals.list",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		HaplotypeCaller \
		-O {output.vcf} \
		-R {params.ref} \
		-I {input.bam}\
		-L {params.interval} \
		--dbsnp {params.dbsnp138} > {log} 2>&1
		"""

rule mergevcfs:
	input:
		vcf=lambda wildcards: expand(
			"{outpath}/03_variants/haplotypecaller/00_initial_call/01_sub_vcf/{sample}.{chr}.g.vcf.gz", 
			outpath=wildcards.outpath, 
			sample=wildcards.sample, 
			chr=CHROMOSOMES
		),
		tbi=lambda wildcards: expand(
			"{outpath}/03_variants/haplotypecaller/00_initial_call/01_sub_vcf/{sample}.{chr}.g.vcf.gz.tbi", 
			outpath=wildcards.outpath, 
			sample=wildcards.sample, 
			chr=CHROMOSOMES
		)
	output:
		vcf_gz="{outpath}/03_variants/haplotypecaller/00_initial_call/02_merge_rawvcf/{sample}.g.vcf.gz",
		vcf_tbi="{outpath}/03_variants/haplotypecaller/00_initial_call/02_merge_rawvcf/{sample}.g.vcf.gz.tbi"
	params:
		gatk=config['gatk_current_using'],
		input_args=lambda wildcards, input: " ".join(f"--INPUT {vcf}" for vcf in input.vcf),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		'''
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		MergeVcfs \
		{params.input_args} \
		--OUTPUT {output.vcf_gz}
		'''

rule sample_name_map:
	input:
		gvcfs=expand("{outpath}/03_variants/haplotypecaller/00_initial_call/02_merge_rawvcf/{sample}.g.vcf.gz", outpath=Outpath, sample=Sample),
		idx=expand("{outpath}/03_variants/haplotypecaller/00_initial_call/02_merge_rawvcf/{sample}.g.vcf.gz.tbi", outpath=Outpath, sample=Sample)
	output:
		"{outpath}/03_variants/haplotypecaller/00_initial_call/02_merge_rawvcf/sample.list.txt"
	log:
		"{outpath}/03_variants/logs/genomics_db_import/sample.list.log"
	params:
		sample_list=lambda wildcards, input: "\n".join([f"{sample}\t{gvcf}" for sample, gvcf in zip(Sample, input.gvcfs)])
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		echo "{params.sample_list}" > {output}
		"""

rule genomics_db_import:
	input:
		"{outpath}/03_variants/haplotypecaller/00_initial_call/02_merge_rawvcf/sample.list.txt"
	output:
		db=directory("{outpath}/03_variants/haplotypecaller/00_initial_call/02_genomicsdb/{chr}/{chr}_gdb")
	log:
		"{outpath}/03_variants/logs/{chr}.genomics_db_import.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval_list=intervals_dir + "/{chr}.intervals.list",
		tmp_dir="{outpath}/03_variants/haplotypecaller/00_initial_call/02_genomicsdb/temp_dir_{chr}",
		gvcf_list=lambda wildcards, input: " ".join([f"-V {gvcf}" for gvcf in input]),
		intervals=lambda wildcards: f"{wildcards.chr}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	shell:
		"""
		mkdir -p {params.tmp_dir}
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		GenomicsDBImport \
		--sample-name-map {input} \
		--genomicsdb-workspace-path {output.db} \
		--tmp-dir {params.tmp_dir} \
		--batch-size 90 \
		-R {params.ref} \
		--max-num-intervals-to-import-in-parallel 25 \
		--intervals {params.intervals} > {log} 2>&1
		"""

# Rule to perform joint genotyping on GenomicsDB workspace
rule genotype_gvcfs:
	input:
		db="{outpath}/03_variants/haplotypecaller/00_initial_call/02_genomicsdb/{chr}/{chr}_gdb"
	output:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/03_joint/{chr}.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/03_joint/{chr}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{chr}.genotype_gvcfs.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		GenotypeGVCFs \
		-R {params.ref} \
		-V gendb://{input.db} \
		-O {output.vcf} \
		--include-non-variant-sites false > {log} 2>&1
		"""

# Rule to merge per-chromosome VCFs into a single cohort VCF
rule merge_vcfs:
	input:
		vcfs=expand("{outpath}/03_variants/haplotypecaller/00_initial_call/03_joint/{chr}.vcf.gz", outpath=Outpath, chr=CHROMOSOMES),
		idx=expand("{outpath}/03_variants/haplotypecaller/00_initial_call/03_joint/{chr}.vcf.gz.tbi", outpath=Outpath, chr=CHROMOSOMES)
	output:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/merge_vcfs.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	params:
		gatk=config['gatk_current_using'],
		vcf_list=lambda wildcards, input: " ".join([f"-I {vcf}" for vcf in input.vcfs]),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		MergeVcfs \
		{params.vcf_list} \
		-O {output.vcf} > {log} 2>&1
		"""

# Rule to perform Variant Quality Score Recalibration (VQSR) for SNPs
rule variant_recalibrator_snp:
	input:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz.tbi"
	output:
		recal="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.snp.recal",
		tranches="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.snp.tranches",
		rscript="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.snp.plots.R"
	log:
		"{outpath}/03_variants/haplotypecaller/00_initial_call/logs/variant_recalibrator_snp.log"
	params:
		ref=config['reference'],
		hapmap=config["hapmap"],
		omni=config["omni"],
		g1000_known_indels=config["g1000_known_indels"],
		dbsnp=config["dbsnp138"],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		VariantRecalibrator \
			-R {params.ref} \
			-V {input.vcf} \
			--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
			--resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} \
			--resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.g1000_known_indels} \
			--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
			-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
			-mode SNP \
			-O {output.recal} \
			--tranches-file {output.tranches} \
			--rscript-file {output.rscript} \
			--dont-run-rscript \
			--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 > {log} 2>&1
		"""

# Rule to perform VQSR for INDELs
rule variant_recalibrator_indel:
	input:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz.tbi"
	output:
		recal="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.indel.recal",
		tranches="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.indel.tranches",
		rscript="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.indel.plots.R"
	log:
		"{outpath}/03_variants/logs/variant_recalibrator_indel.log"
	params:
		ref=config['reference'],
		mills=config["mills_and_1000g"],
		dbsnp=config["dbsnp138"],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		VariantRecalibrator \
			-R {params.ref} \
			-V {input.vcf} \
			--resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
			--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
			-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
			-mode INDEL \
			-O {output.recal} \
			--tranches-file {output.tranches} \
			--rscript-file {output.rscript} \
			--dont-run-rscript \
			--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
			--max-gaussians 4 > {log} 2>&1
		"""

# Rule to apply VQSR for SNPs
rule apply_vqsr_snp:
	input:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/04_merge_vcf/all.vcf.gz.tbi",
		recal="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.snp.recal",
		tranches="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.snp.tranches",
	output:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.snp.recalibrated.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.snp.recalibrated.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/apply_vqsr_snp/all.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		truth_sensitivity=['truth_sensitivity'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		ApplyVQSR \
			-R {params.ref} \
			-V {input.vcf} \
			-O {output.vcf} \
			--recal-file {input.recal} \
			--tranches-file {input.tranches} \
			-mode SNP \
			--truth-sensitivity-filter-level {params.truth_sensitivity} > {log} 2>&1
		"""

# Rule to apply VQSR for INDELs
rule apply_vqsr_indel:
	input:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.snp.recalibrated.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.snp.recalibrated.vcf.gz.tbi",
		recal="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.indel.recal",
		tranches="{outpath}/03_variants/haplotypecaller/00_initial_call/05_vqsr/all.indel.tranches"
	output:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.recalibrated.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.recalibrated.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/apply_vqsr_indel.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		truth_sensitivity=['truth_sensitivity'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		ApplyVQSR \
			-R {params.ref} \
			-V {input.vcf} \
			-O {output.vcf} \
			--recal-file {input.recal} \
			--tranches-file {input.tranches} \
			-mode INDEL \
			--truth-sensitivity-filter-level {params.truth_sensitivity} > {log} 2>&1
		"""

rule split_vcf:
	input:
		vcf="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.recalibrated.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/00_initial_call/06_apply_vqsr/all.recalibrated.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/haplotypecaller/01_raw/{sample}.{chr}.vcf.gz",
		idx="{outpath}/03_variants/haplotypecaller/01_raw/{sample}.{chr}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/hp_split_vcf.{sample}.{chr}.log"
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
