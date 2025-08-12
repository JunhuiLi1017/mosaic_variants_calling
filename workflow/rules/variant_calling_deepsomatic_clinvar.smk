rule deepsomatic_filter:
	input:
		"{outpath}/03_variants/03_deepsomatic/{sample}/{sample}.{chr}.{ref_version}.output.vcf.gz"
	output:
		vcf="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/01_filter/{sample}.{chr}.{ref_version}.pass.output.vcf.gz"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{ref_version}.deepsomatic_filter.log"
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
		vcf="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/01_filter/{sample}.{chr}.{ref_version}.pass.output.vcf.gz"
	output:
		txt="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/02_annovar_clinvar/{sample}.{chr}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{ref_version}.deepsomatic.annotate_clinvar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/02_annovar_clinvar/{sample}.{chr}",
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

rule merge_clinvar:
	input:
		anno=lambda wildcards: expand(
			"{outpath}/03_variants/03_deepsomatic/00_somatic_anno/02_annovar_clinvar/{sample}.{chr}.{ref_version}_multianno.txt", 
			outpath=wildcards.outpath, 
			sample=wildcards.sample, 
			ref_version=wildcards.ref_version,
			chr=chromosomes
		)
	output:
		txt="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/02_annovar_clinvar/{sample}.{ref_version}_multianno.txt"
	params:
		ref_version=config['ref_version'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		cat <(cat {input.anno} | grep "Chr" | sort | uniq) <(cat {input.anno} | grep -v "Chr") > {output.txt}
		"""

rule annotate_rcnv_gnomadlof:
	input:
		tier_anno="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/02_annovar_clinvar/{sample}.{ref_version}_multianno.txt"
	output:
		sub="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/02_annovar_clinvar/score/{sample}.SNV.tier.{ref_version}.exonic_splicing_multianno.txt",
		txt="{outpath}/03_variants/03_deepsomatic/00_somatic_anno/02_annovar_clinvar/score/{sample}.SNV.tier.{ref_version}.rcnv_gnomadlof_multianno.txt"
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
		cat <(awk '{{if($7=="exonic"){{print $0}}}}' {input.tier_anno} | grep -E 'nonsynonymous|stop') <(awk '{{if($7=="splicing"){{print $0}}}}' {input.tier_anno}) > {output.sub}
		cat <(paste <(head -n 1 {input.tier_anno}) <(head -n 1 {params.gnomad_LoF}) <(head -n 1 {params.rCNV_gene_score})) <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$8]){{print $0"\t"c[$8]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.gnomad_LoF} <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$8]){{print $0"\t"c[$8]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.rCNV_gene_score} {output.sub})) > {output.txt}
		"""
