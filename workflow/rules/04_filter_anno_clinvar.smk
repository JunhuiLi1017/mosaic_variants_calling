# Rule 2: Select SNPs
rule split_multiallelic:
	input:
		vcf="{outpath}/03_variants/{caller}/01_raw/{sample}.{chr}.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/01_raw/{sample}.{chr}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/{caller}/02_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/02_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{caller}.split_multiallelic.log"
	params:
		interval=intervals_dir + "/{chr}.intervals.list"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		"""
		bcftools norm \
		  -m - \
		  -o {output.vcf} \
		  -Oz \
		  {input.vcf} > {log} 2>&1
		"""

rule pass_filter:
	input:
		vcf="{outpath}/03_variants/{caller}/02_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/02_split_mulalle/{sample}.{chr}.split_mulalle.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/{caller}/03_pass/{sample}.{chr}.hard_filter.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/03_pass/{sample}.{chr}.hard_filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{caller}.hard_filter.log"
	params:
		interval=intervals_dir + "/{chr}.intervals.list"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		"""
		bcftools filter \
		  -i 'FILTER="PASS" && DP>=10' \
		  -o {output.vcf} \
		  -Oz \
		  {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		"""

# Filter common SNPs using dbSNP, gnomAD, and 1000 Genomes
rule filter_common_snps:
	input:
		vcf="{outpath}/03_variants/{caller}/03_pass/{sample}.{chr}.hard_filter.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/03_pass/{sample}.{chr}.hard_filter.vcf.gz.tbi"
	output:
		vcf_dbsn_filter="{outpath}/03_variants/{caller}/04_commonsnp_filter/{sample}.{chr}.comsnp_dbsnp.vcf.gz",
		vcf_gnomad_filter="{outpath}/03_variants/{caller}/04_commonsnp_filter/{sample}.{chr}.comsnp_gnomad.vcf.gz",
		vcf_pon_filter="{outpath}/03_variants/{caller}/04_commonsnp_filter/{sample}.{chr}.comsnp_filter.vcf.gz",
		vcf_pon_filter_tbi="{outpath}/03_variants/{caller}/04_commonsnp_filter/{sample}.{chr}.comsnp_filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{caller}.comsnp_filter.log"
	params:
		interval=intervals_dir + "/{chr}.intervals.list",
		dbsnp=config["dbsnp138"],
		gnomad_0_01=hg38_sub_gnomad211_exome_genome_dir + "/{chr}.hg38_gnomad211_exome_genome_afpop_gt_0.01_chr.sort.txt",
		pon=config["pon"],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		"""
		bcftools annotate -a {params.dbsnp} -c ID {input.vcf} | bcftools view -e 'ID != "."' | bgzip > {output.vcf_dbsn_filter} && tabix -p vcf {output.vcf_dbsn_filter}
		bcftools view -T ^{params.gnomad_0_01} {output.vcf_dbsn_filter} | bgzip > {output.vcf_gnomad_filter} && tabix -p vcf {output.vcf_gnomad_filter}
		bcftools view -T ^{params.pon} {output.vcf_gnomad_filter} | bgzip > {output.vcf_pon_filter} && tabix -p vcf {output.vcf_pon_filter}
		"""

rule merge_chr:
	input:
		vcf=lambda wildcards: expand("{outpath}/03_variants/{caller}/04_commonsnp_filter/{sample}.{chr}.comsnp_filter.vcf.gz",
		outpath=wildcards.outpath,
		caller=wildcards.caller,
		sample=wildcards.sample,
		chr=CHROMOSOMES),
		tbi=lambda wildcards: [f"{wildcards.outpath}/03_variants/{wildcards.caller}/04_commonsnp_filter/{wildcards.sample}.{chr}.comsnp_filter.vcf.gz.tbi"
		for chr in CHROMOSOMES]
	output:
		vcf="{outpath}/03_variants/{caller}/05_merge/{sample}.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/05_merge/{sample}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{caller}.merge_chr.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		"""
		bcftools concat -a {input.vcf} | bcftools sort | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

rule annotate_clinvar:
	input:
		vcf="{outpath}/03_variants/{caller}/05_merge/{sample}.vcf.gz"
	output:
		txt="{outpath}/03_variants/{caller}/06_annovar/01_annovar/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.{caller}.annovar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/{caller}/06_annovar/01_annovar/{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["terra_perl_anno"]
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

rule annotate_rcnv_gnomadlof_germ:
	input:
		tier_anno="{outpath}/03_variants/{caller}/06_annovar/01_annovar/{sample}.{ref_version}_multianno.txt"
	output:
		txt="{outpath}/03_variants/{caller}/06_annovar/02_score/{sample}.{ref_version}.rcnv_gnomadlof_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.{caller}.rcnv_lof.log"
	params:
		ref_version=config['ref_version'],
		gnomad_LoF=config['gnomad_LoF'],
		rCNV_gene_score=config['rCNV_gene_score'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		sub="{outpath}/03_variants/{caller}/06_annovar/02_score/{sample}.{ref_version}.exonic_splicing_multianno.txt"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		"""
		cat <(awk '{{if($7=="exonic"){{print $0}}}}' {input.tier_anno} | grep -E 'nonsynonymous|stop') <(awk '{{if($7=="splicing"){{print $0}}}}' {input.tier_anno}) > {params.sub}
		cat <(paste <(head -n 1 {input.tier_anno}) <(head -n 1 {params.gnomad_LoF}) <(head -n 1 {params.rCNV_gene_score})) <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$8]){{print $0"\t"c[$8]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.gnomad_LoF} <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$8]){{print $0"\t"c[$8]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.rCNV_gene_score} {params.sub})) > {output.txt}
		rm {params.sub}
		"""

rule reformat_rcnv_gnomadlof_germ:
	input:
		txt="{outpath}/03_variants/{caller}/06_annovar/02_score/{sample}.{ref_version}.rcnv_gnomadlof_multianno.txt"
	output:
		txt="{outpath}/03_variants/{caller}/06_annovar/03_reformat_anno/{sample}.{ref_version}.rcnv_gnomadlof_multianno.reformat.txt"
	log:
		"{outpath}/03_variants/{sample}.{caller}.{ref_version}.m2.reformat.log"
	params:
		reformat_script=config['reformat_script'],
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["python_alpine3.21"]
	shell:
		"""
		python {params.reformat_script} {input.txt} {output.txt} --info-column "Otherinfo12" --gt-column "Otherinfo13" > {log} 2>&1
		"""

rule reprot_germ:
	input:
		expand("{outpath}/03_variants/{caller}/06_annovar/03_reformat_anno/{sample}.{ref_version}.rcnv_gnomadlof_multianno.reformat.txt",
			outpath=Outpath,
			caller=Callers,
			sample=Sample,
			ref_version=Ref_version)
	output:
		"{outpath}/03_variants/variant_report.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		wc -l {input} > {output}
		'''