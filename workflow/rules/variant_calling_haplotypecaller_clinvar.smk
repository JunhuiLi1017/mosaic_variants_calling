rule split_multiallelic:
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/00_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/00_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.hc.split_multiallelic.log"
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


# Rule 2: Select SNPs
rule select_snps:
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/00_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/00_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_snp/{sample}.{chr}.snp.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_snp/{sample}.{chr}.snp.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.select_snps.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval=intervals_dir + "/{chr}.intervals.list",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
		SelectVariants \
		-O {output.vcf} \
		-R {params.ref} \
		-V {input.vcf}	\
		--select-type SNP > {log} 2>&1
		conda deactivate
		"""

# Rule 3: Select Indels
rule select_indels:
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/sub_vcf/{sample}/{sample}.{chr}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_indel/{sample}.{chr}.indel.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_indel/{sample}.{chr}.indel.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.select_indels.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval=intervals_dir + "/{chr}.intervals.list",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
		SelectVariants \
		-O {output.vcf} \
		-R {params.ref} \
		-V {input.vcf}	\
		--select-type INDEL > {log} 2>&1
		conda deactivate
		"""

# Rule 4: Hard-filter SNPs
rule filter_snps:
	message:
		"GATK: hard filtering germline short variants: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants"
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_snp/{sample}.{chr}.snp.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_snp/{sample}.{chr}.snp.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/02_sub_vcf_snp_hard_filter/{sample}.{chr}.snp.filter.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/02_sub_vcf_snp_hard_filter/{sample}.{chr}.snp.filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.filter_snps.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval=intervals_dir + "/{chr}.intervals.list",
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
		bcftools filter \
		-e 'QD < 2.0 || QUAL < 30.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -3.0 || INFO/DP < 10' \
		{input.vcf} \
		-s "HARD_FILTER_SNP" \
		--threads {threads} | \
		bcftools filter -i 'FILTER="PASS"' | \
		bcftools sort | \
		bgzip -c > {output.vcf} 2> {log}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

# Rule 5: Hard-filter INDELs
rule filter_indels:
	message:
		"GATK: hard filtering germline short variants: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants"
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_indel/{sample}.{chr}.indel.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/01_sub_vcf_indel/{sample}.{chr}.indel.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/02_sub_vcf_indel_hard_filter/{sample}.{chr}.indel.filter.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/02_sub_vcf_indel_hard_filter/{sample}.{chr}.indel.filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.filter_indels.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval=intervals_dir + "/{chr}.intervals.list",
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
		bcftools filter \
		-e 'QD < 2.0 || QUAL < 30.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || INFO/DP < 10' \
		{input.vcf} \
		-s "HARD_FILTER_INDEL" \
		--threads {threads} | \
		bcftools filter -i 'FILTER="PASS"' | \
		bcftools sort | \
		bgzip -c > {output.vcf} 2> {log}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

# Rule 6: Merge filtered SNPs and Indels
rule merge_variants:
	input:
		snps="{outpath}/03_variants/02_haplotypecaller/02_sub_vcf_snp_hard_filter/{sample}.{chr}.snp.filter.vcf.gz",
		indels="{outpath}/03_variants/02_haplotypecaller/02_sub_vcf_indel_hard_filter/{sample}.{chr}.indel.filter.vcf.gz"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/03_sub_vcf_variant_hard_filter/{sample}.{chr}.hard_filter.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/03_sub_vcf_variant_hard_filter/{sample}.{chr}.hard_filter.vcf.gz.tbi"
	params:
		interval=intervals_dir + "/{chr}.intervals.list"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools concat -a {input.snps} {input.indels} | bcftools sort | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

# Rule 7: Filter common SNPs using dbSNP, gnomAD, and 1000 Genomes
rule filter_common_dbsnp:
	message:
		"""
		bcftools view -e 'ID != "." || AF>0.01' gnomad.vcf | \ ## this will automatically remove variants with af > 0.01 but not present in gnomda.
		"""
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/03_sub_vcf_variant_hard_filter/{sample}.{chr}.hard_filter.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/03_sub_vcf_variant_hard_filter/{sample}.{chr}.hard_filter.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_dbsnp.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_dbsnp.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.filter_indels_dbsnp.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		gnomad_0_001=hg38_sub_gnomad211_exome_genome_dir + "/{chr}.hg38_gnomad211_exome_genome_afpop_gt_0.01_chr.sort.txt",
		dbsnp=config["dbsnp138"],
		pon=config["pon"],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools	
		bcftools view -e 'ID != "."' {input.vcf} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

rule filter_common_pon:
	message:
		"""
		bcftools view -e 'ID != "." || AF>0.01' gnomad.vcf | \ ## this will automatically remove variants with af > 0.01 but not present in gnomda.
		"""
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_dbsnp.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_dbsnp.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_pon.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_pon.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.filter_indels_pon.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		gnomad_0_001=hg38_sub_gnomad211_exome_genome_dir + "/{chr}.hg38_gnomad211_exome_genome_afpop_gt_0.01_chr.sort.txt",
		dbsnp=config["dbsnp138"],
		pon=config["pon"],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools	
		bcftools view -T ^{params.pon} {input.vcf} |  bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

rule filter_common_gnomad:
	message:
		"""
		bcftools view -e 'ID != "." || AF>0.01' gnomad.vcf | \ ## this will automatically remove variants with af > 0.01 but not present in gnomda.
		"""
	input:
		vcf="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_pon.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_pon.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_gnomad.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.{chr}.common_filter_gnomad.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.filter_indels_gnomad.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		gnomad_0_001=hg38_sub_gnomad211_exome_genome_dir + "/{chr}.hg38_gnomad211_exome_genome_afpop_gt_0.01_chr.sort.txt",
		dbsnp=config["dbsnp138"],
		pon=config["pon"],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools	
		bcftools view -T ^{params.gnomad_0_001} {input.vcf} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""
		
rule merge_common_filter:
	input:
		vcf=expand("{outpath}/03_variants/02_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{{sample}}.{chr}.common_filter_gnomad.vcf.gz", outpath=outpath, chr=chromosomes)
	output:
		vcf="{outpath}/03_variants/02_haplotypecaller/05_merge_common_filter/{sample}.common_filter.vcf.gz",
		tbi="{outpath}/03_variants/02_haplotypecaller/05_merge_common_filter/{sample}.common_filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.merge_common_filter.log"
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
		vcf="{outpath}/03_variants/02_haplotypecaller/05_merge_common_filter/{sample}.common_filter.vcf.gz"
	output:
		txt="{outpath}/03_variants/02_haplotypecaller/06_annovar/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.06_annovar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/02_haplotypecaller/06_annovar/{sample}",
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
		-protocol refGene,dbnsfp42a,clinvar_20240917 \
		-operation g,f,f \
		-nastring . \
		-vcfinput
		"""

rule annotate_rcnv_gnomadlof:
	input:
		tier_anno="{outpath}/03_variants/02_haplotypecaller/06_annovar/{sample}.{ref_version}_multianno.txt"
	output:
		sub="{outpath}/03_variants/02_haplotypecaller/06_annovar/score/{sample}.SNV.tier.{ref_version}.exonic_splicing_multianno.txt",
		txt="{outpath}/03_variants/02_haplotypecaller/06_annovar/score/{sample}.SNV.tier.{ref_version}.rcnv_gnomadlof_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.m2.rcnv_lof.log"
	params:
		ref_version=config['ref_version'],
		gnomad_LoF=config['gnomad_LoF'],
		rCNV_gene_score=config['rCNV_gene_score'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		"""
		cat <(awk '{{if($6=="exonic"){{print $0}}}}' {input.tier_anno} | grep -E 'nonsynonymous|stop') <(awk '{{if($6=="splicing"){{print $0}}}}' {input.tier_anno}) > {output.sub}
		cat <(paste <(head -n 1 {input.tier_anno}) <(head -n 1 {params.gnomad_LoF}) <(head -n 1 {params.rCNV_gene_score})) <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$7]){{print $0"\t"c[$7]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.gnomad_LoF} <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$7]){{print $0"\t"c[$7]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.rCNV_gene_score} {output.sub})) > {output.txt}
		"""
