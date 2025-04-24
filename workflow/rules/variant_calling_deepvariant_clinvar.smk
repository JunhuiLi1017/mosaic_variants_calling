# Rule 2: Select SNPs
rule split_multiallelic:
	input:
		vcf="{outpath}/03_variants/04_deepvariant/{sample}/{sample}.{chr}.{ref_version}.output.vcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/{sample}/{sample}.{chr}.{ref_version}.output.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/04_deepvariant/01_split_mulalle/{sample}.{chr}.{ref_version}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/01_split_mulalle/{sample}.{chr}.{ref_version}.split_mulalle.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{ref_version}.split_multiallelic.log"
	params:
		interval=intervals_dir + "/{chr}.intervals.list"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools norm \
		  -m - \
		  -o {output.vcf} \
		  -Oz \
		  {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		"""

rule hard_filter:
	input:
		vcf="{outpath}/03_variants/04_deepvariant/01_split_mulalle/{sample}.{chr}.{ref_version}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/01_split_mulalle/{sample}.{chr}.{ref_version}.split_mulalle.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/04_deepvariant/02_hard_filter/{sample}.{chr}.{ref_version}.hard_filter.vcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/02_hard_filter/{sample}.{chr}.{ref_version}.hard_filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{ref_version}.hard_filter.log"
	params:
		interval=intervals_dir + "/{chr}.intervals.list"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools filter \
		  -i 'FILTER="PASS" && DP>=10' \
		  -o {output.vcf} \
		  -Oz \
		  {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

# Filter common SNPs using dbSNP, gnomAD, and 1000 Genomes
rule filter_common_snps:
	input:
		vcf="{outpath}/03_variants/04_deepvariant/02_hard_filter/{sample}.{chr}.{ref_version}.hard_filter.vcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/02_hard_filter/{sample}.{chr}.{ref_version}.hard_filter.vcf.gz.tbi"
	output:
		vcf_dbsn_filter="{outpath}/03_variants/04_deepvariant/03_commonsnp_filter/{sample}.{chr}.{ref_version}.comsnp_dbsnp.vcf.gz",
		vcf_gnomad_filter="{outpath}/03_variants/04_deepvariant/03_commonsnp_filter/{sample}.{chr}.{ref_version}.comsnp_gnomad.vcf.gz",
		vcf_pon_filter="{outpath}/03_variants/04_deepvariant/03_commonsnp_filter/{sample}.{chr}.{ref_version}.comsnp_filter.vcf.gz",
		vcf_pon_filter_tbi="{outpath}/03_variants/04_deepvariant/03_commonsnp_filter/{sample}.{chr}.{ref_version}.comsnp_filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{ref_version}.comsnp_filter.log"
	params:
		interval=intervals_dir + "/{chr}.intervals.list",
		dbsnp=config["dpsnp138_gz"],
		gnomad_0_01=hg38_sub_gnomad211_exome_genome_dir + "/{chr}.hg38_gnomad211_exome_genome_afpop_gt_0.01_chr.sort.txt",
		pon=config["pon"],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools annotate -a {params.dbsnp} -c ID {input.vcf} | bcftools view -e 'ID != "."' | bgzip > {output.vcf_dbsn_filter} && tabix -p vcf {output.vcf_dbsn_filter}
		bcftools view -T ^{params.gnomad_0_01} {output.vcf_dbsn_filter} | bgzip > {output.vcf_gnomad_filter} && tabix -p vcf {output.vcf_gnomad_filter}
		bcftools view -T ^{params.pon} {output.vcf_gnomad_filter} | bgzip > {output.vcf_pon_filter} && tabix -p vcf {output.vcf_pon_filter}
		conda deactivate
		"""

rule merge_common_filter:
	input:
		vcf=lambda wildcards: [f"{wildcards.outpath}/03_variants/04_deepvariant/03_commonsnp_filter/{wildcards.sample}.{chr}.{wildcards.ref_version}.comsnp_filter.vcf.gz" for chr in chromosomes],
		tbi=lambda wildcards: [f"{wildcards.outpath}/03_variants/04_deepvariant/03_commonsnp_filter/{wildcards.sample}.{chr}.{wildcards.ref_version}.comsnp_filter.vcf.gz.tbi" for chr in chromosomes]
	output:
		vcf="{outpath}/03_variants/04_deepvariant/05_merge_common_filter/{sample}.common_filter.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/05_merge_common_filter/{sample}.common_filter.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.merge_common_filter.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools concat -a {input.vcf} | bcftools sort | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

rule annotate_clinvar:
	input:
		vcf="{outpath}/03_variants/04_deepvariant/05_merge_common_filter/{sample}.common_filter.{ref_version}.vcf.gz"
	output:
		txt="{outpath}/03_variants/04_deepvariant/06_annovar/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.06_annovar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/04_deepvariant/06_annovar/{sample}",
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