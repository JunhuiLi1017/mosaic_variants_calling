# Merged Snakemake file for variant calling with Mutect2 and MosaicForecast
# Supports both Refine and Phase models

rule Mutect2:
	input:
		"{outpath}/02_map/08_bqsr/{sample}.bam"
	output:
		raw_vcf="{outpath}/03_variants/01_mutect2/01_somatic/01_result_sub/{sample}.{individual_chr}.mt2pon.vcf.gz",
		raw_stat="{outpath}/03_variants/01_mutect2/01_somatic/01_result_sub/{sample}.{individual_chr}.mt2pon.vcf.gz.stats"
	log:
		"{outpath}/03_variants/logs/{sample}.{individual_chr}.mutect.log"
	params:
		pon=config['pon'],
		af_only_gnomad=config['af_only_gnomad'],
		ref=config['reference'],
		sample="{sample}",
		interval_list=intervals_dir + "/{individual_chr}.intervals.list",
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		'''
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} Mutect2 \
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
		lambda wildcards: [f"{wildcards.outpath}/03_variants/01_mutect2/01_somatic/01_result_sub/{wildcards.sample}.{chr}.mt2pon.vcf.gz" for chr in CHROMOSOMES]
	output:
		args="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.args",
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.merged.vcf.gz",
		idx="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.merged.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.merge.logs"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		ls {input} > {output.args}
		vcf-concat -f {output.args}  |  bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		'''

rule merge_mutectstats:
	input:
		lambda wildcards: [f"{wildcards.outpath}/03_variants/01_mutect2/01_somatic/01_result_sub/{wildcards.sample}.{chr}.mt2pon.vcf.gz.stats" for chr in CHROMOSOMES]
	output:
		"{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.merged.vcf.gz.stats"
	params:
		list_para = lambda wildcards: ' '.join([f"--stats {wildcards.outpath}/03_variants/01_mutect2/01_somatic/01_result_sub/{wildcards.sample}.{chr}.mt2pon.vcf.gz.stats" for chr in CHROMOSOMES]),
		gatk=config['gatk_current_using']
	log:
		"{outpath}/03_variants/logs/{sample}.MergeMutectStats.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		'''
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
		MergeMutectStats \
		{params.list_para} \
		-O {output} > {log} 2>&1
		'''

rule FilterMutectCall:
	input:
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.merged.vcf.gz",
		stats="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.merged.vcf.gz.stats",
		index="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.merged.vcf.gz.tbi"
	output:
		"{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz"
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
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		'''
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
		FilterMutectCalls -R {params.ref} -V {input.vcf} -O {output} > {log} 2>&1
		'''

##--------------------------
## somatic filter and annotation
##--------------------------

rule mutect2_pass:
	input:
		"{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz"
	output:
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/02_pass/{sample}.{ref_version}.pass.output.vcf.gz"
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
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/02_pass/{sample}.{ref_version}.pass.output.vcf.gz"
	output:
		txt="{outpath}/03_variants/01_mutect2/01_somatic/03_annotation/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.m2.annotate_clinvar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/01_mutect2/01_somatic/03_annotation/{sample}",
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

##--------------------------
## mosaicforecast - unified for both Refine and Phase models
##--------------------------

# Helper function to determine output paths based on model
def get_mosaicforecast_paths(wildcards):
    if wildcards.model == "Phase":
        base_path = f"{wildcards.outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast"
    else:  # Refine model
        base_path = f"{wildcards.outpath}/03_variants/01_mutect2/02_mosaicforecast"
    return base_path

rule MT2_initial_filter:
	input:
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/02_init_filter/{sample}.mt2pon.AF0.02.bed"
	params:
		sample="{sample}",
		af_threshold= lambda wildcards: "0.03" if units[units['sample'] == wildcards.sample]['pcr_based'].any() else "0.02"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat <(zcat {input.vcf} \
		| grep -v "^#" \
		| grep PASS \
		| gawk '{{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){{print $0}}}}' \
		| cut -f1,2,4,5,10 \
		| sed "s/:/\\t/g" \
		| sed "s/,/\\t/g" \
		| awk '$8>=0.03 && $8<0.4' \
		| grep -v '0|1' \
		| grep -v '1|0') \
		<(zcat {input.vcf} \
		| grep -v "^#" \
		| grep PASS \
		| gawk '{{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){{print $0}}}}' \
		| cut -f1,2,4,5,10 \
		| sed "s/:/\\t/g" \
		| sed "s/,/\\t/g" \
		| awk '$8>={params.af_threshold} && $8<0.4 ' \
		| grep -E "0\|1|1\|0") \
		| cut -f 1-4,6-8 \
		| awk '{{OFS="\\t";print $1,$2-1,$2,$3,$4,\"{params.sample}\",$5,$6,$7}}' > {output}
		'''

rule repeat_filter:
	input:
		mt2pon="{outpath}/03_variants/01_mutect2/02_mosaicforecast/02_init_filter/{sample}.mt2pon.AF0.02.bed"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.noSegDup.bed"
	params:
		segdup=config['SegDup_and_clustered']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	conda:
		"../envs/bedtools.yaml"
	shell:
		'''
		subtractBed -a {input.mt2pon} -b {params.segdup} > {output}
		'''

rule annovar_formatter:
	input:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.noSegDup.bed"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat {input} | awk '{{OFS="\\t\";len=length($4)-length($5);if(len<=0){{print $1,$3,$3,$4,$5,$6}}if(len>0){{print $1,$3,$3+len,$4,$5,$6}}}}'> {output}
		sed -i 's/chr//g' {output}
		'''

rule MAF0_extraction_SNV:
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.SNV.chr.bed"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)==1 && length($5)==1'\'' |awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"

rule MAF0_extraction_INS:
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.INS.chr.bed"
	params:
		allrepeats_forindel=config['allrepeats_forindel']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	conda:
		"../envs/bedtools.yaml"
	shell:
		"subtractBed -a <(cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' |awk '\''length($4)< length($5)'\'') -b {params.allrepeats_forindel} | awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"

rule MAF0_extraction_DEL:
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.DEL.chr.bed"
	params:
		allrepeats_forindel=config['allrepeats_forindel']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	conda:
		"../envs/bedtools.yaml"
	shell:
		"subtractBed -a <(cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)> length($5)'\'') -b {params.allrepeats_forindel} | awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"

rule split_bed:
	input:
		bed_snv=f"{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.SNV.chr.bed",
		bed_ins=f"{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.INS.chr.bed",
		bed_del=f"{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.DEL.chr.bed"
	output:
		bed_chunks_snv="{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.snv.bed",
		bed_chunks_ins="{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.ins.bed",
		bed_chunks_del="{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.del.bed"
	threads:
		resource['resource']['low']['threads']
	params:
		individual_chr="{individual_chr}"
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		# Extract SNVs for this chromosome
		awk -v chr="{params.individual_chr}" '$1=="chr"chr' {input.bed_snv} > {output.bed_chunks_snv}

		# Extract INS for this chromosome 
		awk -v chr="{params.individual_chr}" '$1=="chr"chr' {input.bed_ins} > {output.bed_chunks_ins}

		# Extract DEL for this chromosome
		awk -v chr="{params.individual_chr}" '$1=="chr"chr' {input.bed_del} > {output.bed_chunks_del}
		"""

rule feature_extraction_SNV_chunk:
	input:
		bed="{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.snv.bed",
		bam="{outpath}/02_map/05_bqsr/{sample}.bam",
		bai="{outpath}/02_map/05_bqsr/{sample}.bai"
	output:
		feature="{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/05_feature_each_chromosome/{sample}.SNV.{individual_chr}.features"
	params:
		outdir="{outpath}/02_map/05_bqsr",
		umap=config['umap'],
		ReadLevel_Features_extraction=config['ReadLevel_Features_extraction'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/{sample}.SNV.{individual_chr}.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate py3.7.1
		export PYTHONWARNINGS="ignore" 
		python {params.ReadLevel_Features_extraction} \
		{input.bed} \
		{output.feature} \
		{params.outdir} \
		{params.ref} \
		{params.umap} \
		1 \
		bam > {log} 2>&1 
		conda deactivate
		'''

rule merge_features_SNV:
	input:
		lambda wildcards: [f"{wildcards.outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/05_feature_each_chromosome/{wildcards.sample}.SNV.{chr}.features" for chr in CHROMOSOMES]
	output:
		feature="{outpath}/03_variants/01_mutect2_pon/02_mosaicforecast/05_feature/{sample}.SNV.features"
	threads: 1
	resources:
		mem_mb=1000
	shell:
		'''
		# Get header from first file
		head -n 1 {input[0]} > {output.feature}
		# Append all data lines from all files
		for f in {input}; do
			tail -n +2 "$f" >> {output.feature}
		done
		'''

rule feature_extraction_INS:
	message:
		"since snakemake cannot create a conda env from ../envs/py3.7.1.yaml, but we would like to use --use-conda, so we use conda activate"
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.INS.chr.bed",
		file2="{outpath}/02_map/05_bqsr/{sample}.bam"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/05_feature/{sample}.INS.features"
	params:
		outdir="{outpath}/02_map/05_bqsr/",
		umap=config['umap'],
		ReadLevel_Features_extraction=config['ReadLevel_Features_extraction'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/{sample}.INS.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate py3.7.1
		export PYTHONWARNINGS="ignore" 
		python {params.ReadLevel_Features_extraction} {input.file1} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1
		conda deactivate
		'''

rule feature_extraction_DEL:
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.DEL.chr.bed",
		file2="{outpath}/02_map/05_bqsr/{sample}.bam"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/05_feature/{sample}.DEL.features"
	params:
		outdir="{outpath}/02_map/05_bqsr/",
		umap=config['umap'],
		ReadLevel_Features_extraction=config['ReadLevel_Features_extraction'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/{sample}.DEL.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate py3.7.1
		export PYTHONWARNINGS="ignore" 
		python {params.ReadLevel_Features_extraction} \
		{input.file1} \
		{output} \
		{params.outdir} \
		{params.ref} \
		{params.umap} \
		1 \
		bam > {log} 2>&1
		conda deactivate
		'''

rule g1000_avail_acess_filter_SNV:
	input:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/05_feature/{sample}.SNV.features"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/06_1000g_avail_access/{sample}.SNV.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		cat <(head -n 1 {input}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input} | cut -f 2-3) {input}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_INS:
	input:
		feature="{outpath}/03_variants/01_mutect2/02_mosaicforecast/05_feature/{sample}.INS.features"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/06_1000g_avail_access/{sample}.INS.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		cat <(head -n 1 {input.feature}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_DEL:
	input:
		feature="{outpath}/03_variants/01_mutect2/02_mosaicforecast/05_feature/{sample}.DEL.features"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/06_1000g_avail_access/{sample}.DEL.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat <(head -n 1 {input.feature}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		''' 

rule Prediction_SNV:
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/06_1000g_avail_access/{sample}.SNV.no1000g.features"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.predictions"
	params:
		refine_beta=config['refine_beta'],
		phase_model=config['phase_model'],
		model="{model}",
		prediction_r_script=config['prediction_r_script']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		if [ "{params.model}" = "Phase" ]; then
			Rscript {params.prediction_r_script} {input.file1} {params.phase_model} Phase {output}
		else
			Rscript {params.prediction_r_script} {input.file1} {params.refine_beta} {params.model} {output}
		fi
		'''

rule Prediction_INS:
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/06_1000g_avail_access/{sample}.INS.no1000g.features"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.predictions"
	params:
		model="{model}",
		refine_beta=config['prediction_indel_model'],
		phase_model=config['phase_model'],
		prediction_r_script=config['prediction_r_script']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		if [ "{params.model}" = "Phase" ]; then
			Rscript {params.prediction_r_script} {input.file1} {params.phase_model} Phase {output}
		else
			Rscript {params.prediction_r_script} {input.file1} {params.refine_beta} {params.model} {output}
		fi
		'''

rule Prediction_DEL:
	input:
		file1="{outpath}/03_variants/01_mutect2/02_mosaicforecast/06_1000g_avail_access/{sample}.DEL.no1000g.features"
	output:
		"{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.predictions"
	params:
		model="{model}",
		refine_beta=config['prediction_indel_model'],
		phase_model=config['phase_model'],
		prediction_r_script=config['prediction_r_script']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		if [ "{params.model}" = "Phase" ]; then
			Rscript {params.prediction_r_script} {input.file1} {params.phase_model} Phase {output}
		else
			Rscript {params.prediction_r_script} {input.file1} {params.refine_beta} {params.model} {output}
		fi
		'''

rule extract_bed_snv:
	input:
		prediction_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.predictions"
	output:
		bed_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.{ref_version}.bed.temp_file",
		model="{model}"
	shell:
		'''
		if [ "{params.model}" = "Phase" ]; then
			awk '$35=="hap=3"{{print $0}}' {input.prediction_snv} > {params.temp_file}
		else
			awk '$35=="mosaic"{{print $0}}' {input.prediction_snv} > {params.temp_file}
		fi
		if [ -s {params.temp_file} ]; then
			cat {params.temp_file} | cut -f 1 | cut -d"~" -f2-5 | sed 's/~/\t/g' | grep -v "id" | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n > {output.bed_snv}
		else
			: > {output.bed_snv}
		fi
		rm {params.temp_file}
		'''

rule extract_bed_ins:
	input:
		prediction_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.predictions"
	output:
		bed_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.{ref_version}.bed.temp_file",
		model="{model}"
	shell:
		'''
		if [ "{params.model}" = "Phase" ]; then
			awk '$35=="hap=3"{{print $0}}' {input.prediction_ins} > {params.temp_file}
		else
			awk '$35=="hap=3"{{print $0}}' {input.prediction_ins} > {params.temp_file}
		fi
		if [ -s {params.temp_file} ]; then
			cat {params.temp_file} | cut -f 1 | cut -d"~" -f2-5 | sed 's/~/\t/g' | grep -v "id" | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n > {output.bed_ins}
		else
			: > {output.bed_ins}
		fi
		rm {params.temp_file}
		'''

rule extract_bed_del:
	input:
		prediction_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.predictions"
	output:
		bed_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.{ref_version}.bed.temp_file",
		model="{model}"
	shell:
		'''
		if [ "{params.model}" = "Phase" ]; then
			awk '$35=="hap=3"{{print $0}}' {input.prediction_del} > {params.temp_file}
		else
			awk '$35=="hap=3"{{print $0}}' {input.prediction_del} > {params.temp_file}
		fi
		if [ -s {params.temp_file} ]; then
			cat {params.temp_file} | cut -f 1 | cut -d"~" -f2-5 | sed 's/~/\t/g' | grep -v "id" | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n > {output.bed_del}
		else
			: > {output.bed_del}
		fi
		rm {params.temp_file}
		'''

rule extrac_SNVsubvcf:
	input:
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz",
		bed="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.{ref_version}.bed"
	output:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
		dpaf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.{ref_version}.Geno.DP.AF.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		# Check if BED file is empty
		if [ -s {input.bed} ]; then
			# Non-empty BED: filter VCF and generate outputs
			bcftools view -T {input.bed} {input.vcf} -o {output.vcf} -Oz
			zcat {output.vcf} | grep -v "#" | cut -f 1,2,4,5,10 | sed "s/:/\t/g" | \
				cut -f 1-8 | sed 's/,/\t/g' | \
				awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		else
			# Empty BED: create empty VCF with header and empty DP/AF table
			echo "Warning: No mosaic variants found, creating empty VCF and DP/AF files" >&2
			bcftools view -h {input.vcf} | bgzip > {output.vcf}
			echo -e "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT" > {output.dpaf}
		fi
		tabix -p vcf {output.vcf}
		conda deactivate		
		'''

rule extrac_INSsubvcf:
	input:
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz",
		bed="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.{ref_version}.bed"
	output:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.{ref_version}.vcf.gz",
		dpaf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.{ref_version}.Geno.DP.AF.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		# Check if BED file is empty
		if [ -s {input.bed} ]; then
			# Non-empty BED: filter VCF and generate outputs
			bcftools view -T {input.bed} {input.vcf} -o {output.vcf} -Oz
			zcat {output.vcf} | grep -v "#" | cut -f 1,2,4,5,10 | sed "s/:/\t/g" | \
				cut -f 1-8 | sed 's/,/\t/g' | \
				awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		else
			# Empty BED: create empty VCF with header and empty DP/AF table
			echo "Warning: No mosaic variants found, creating empty VCF and DP/AF files" >&2
			bcftools view -h {input.vcf} | bgzip > {output.vcf}
			echo -e "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT" > {output.dpaf}
		fi
		tabix -p vcf {output.vcf}
		conda deactivate
		'''

rule extrac_DELsubvcf:
	input:
		vcf="{outpath}/03_variants/01_mutect2/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz",
		bed="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.{ref_version}.bed"
	output:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
		dpaf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.{ref_version}.Geno.DP.AF.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	log:
		"{outpath}/logs/extract_DELsubvcf/{sample}.{cov}.{ref_version}.log"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		# Check if BED file is empty
		if [ ! -s "{input.bed}" ]; then
			echo "No variants: creating empty outputs" >> {log}
			bcftools view -h "{input.vcf}" | bgzip > "{output.vcf}" 2>> {log} || true
			echo -e "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT" > "{output.dpaf}" 2>> {log}
		else
			# Process non-empty BED
			bcftools view -T "{input.bed}" "{input.vcf}" -o "{output.vcf}" -Oz 2>> {log}
			zcat "{output.vcf}" | grep -v "#" | cut -f 1,2,4,5,10 | sed "s/:/\t/g" | \
				cut -f 1-8 | sed 's/,/\t/g' | \
				awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > "{output.dpaf}" 2>> {log}
		fi
		tabix -p vcf {output.vcf}
		conda deactivate
		''' 

rule anno_gnomAD_dbsbp_snv:
	input:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
	output:
		vcf_anno_dbsnp="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.SNV.anno.dbsnp.{ref_version}.vcf.gz",
		vcf_anno_exome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.SNV.anno.exome.{ref_version}.vcf.gz",
		vcf_anno_genome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_anno_genome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{cov}.SNV.{ref_version}.dbsnp.gnomad.logs"
	params:
		dbsnp=config["dbsnp138"],
		gnomad_exome_header=config['gnomad_exome_header'],
		gnomad_genome_header=config['gnomad_genome_header'],
		gnomad_genome=config['gnomad_genome'],
		gnomad_exome=config['gnomad_exome']	
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools annotate \
			-a {params.dbsnp} \
			-c ID \
			{input.vcf} \
			-o {output.vcf_anno_dbsnp} \
			-Oz > {log} 2>&1
		tabix -p vcf {output.vcf_anno_dbsnp}

		# Annotate with exome data
		bcftools annotate \
			-a {params.gnomad_exome} \
			-h {params.gnomad_exome_header} \
			-c CHROM,POS,POS,REF,ALT,AF_exome,AF_popmax_exome \
			{output.vcf_anno_dbsnp} \
			-o {output.vcf_anno_exome} \
			-Oz > {log} 2>&1
		tabix -p vcf {output.vcf_anno_exome}

		# Annotate with genome data
		bcftools annotate \
			-a {params.gnomad_genome} \
			-h {params.gnomad_genome_header} \
			-c CHROM,POS,POS,REF,ALT,AF_genome,AF_popmax_genome \
			{output.vcf_anno_exome} \
			-o {output.vcf_anno_genome} \
			-Oz >> {log} 2>&1
		tabix -p vcf {output.vcf_anno_genome}
		'''

rule anno_gnomAD_dbsbp_ins:
	input:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.INS.{ref_version}.vcf.gz",
	output:
		vcf_anno_dbsnp="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.INS.anno.dbsnp.{ref_version}.vcf.gz",
		vcf_anno_exome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.INS.anno.exome.{ref_version}.vcf.gz",
		vcf_anno_genome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_anno_genome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{cov}.INS.{ref_version}.dbsnp.gnomad.logs"
	params:
		dbsnp=config["dbsnp138"],
		gnomad_exome_header=config['gnomad_exome_header'],
		gnomad_genome_header=config['gnomad_genome_header'],
		gnomad_genome=config['gnomad_genome'],
		gnomad_exome=config['gnomad_exome']	
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools annotate \
			-a {params.dbsnp} \
			-c ID \
			{input.vcf} \
			-o {output.vcf_anno_dbsnp} \
			-O z > {log} 2>&1
		tabix -p vcf {output.vcf_anno_dbsnp}

		# Annotate with exome data
		bcftools annotate \
			-a {params.gnomad_exome} \
			-h {params.gnomad_exome_header} \
			-c CHROM,POS,POS,REF,ALT,AF_exome,AF_popmax_exome \
			{output.vcf_anno_dbsnp} \
			-o {output.vcf_anno_exome} \
			-O z > {log} 2>&1
		tabix -p vcf {output.vcf_anno_exome}

		# Annotate with genome data
		bcftools annotate \
			-a {params.gnomad_genome} \
			-h {params.gnomad_genome_header} \
			-c CHROM,POS,POS,REF,ALT,AF_genome,AF_popmax_genome \
			{output.vcf_anno_exome} \
			-o {output.vcf_anno_genome} \
			-O z >> {log} 2>&1
		tabix -p vcf {output.vcf_anno_genome}
		'''

rule anno_gnomAD_dbsbp_del:
	input:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/07_mosaic_prediction/{model}/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
	output:
		vcf_anno_dbsnp="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.DEL.anno.dbsnp.{ref_version}.vcf.gz",
		vcf_anno_exome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.DEL.anno.exome.{ref_version}.vcf.gz",
		vcf_anno_genome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_anno_genome="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{cov}.DEL.{ref_version}.dbsnp.gnomad.logs"
	params:
		dbsnp=config["dbsnp138"],
		gnomad_exome_header=config['gnomad_exome_header'],
		gnomad_genome_header=config['gnomad_genome_header'],
		gnomad_genome=config['gnomad_genome'],
		gnomad_exome=config['gnomad_exome']	
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools annotate \
			-a {params.dbsnp} \
			-c ID \
			{input.vcf} \
			-o {output.vcf_anno_dbsnp} \
			-Oz > {log} 2>&1
		tabix -p vcf {output.vcf_anno_dbsnp}

		# Annotate with exome data
		bcftools annotate \
			-a {params.gnomad_exome} \
			-h {params.gnomad_exome_header} \
			-c CHROM,POS,POS,REF,ALT,AF_exome,AF_popmax_exome \
			{output.vcf_anno_dbsnp} \
			-o {output.vcf_anno_exome} \
			-Oz > {log} 2>&1
		tabix -p vcf {output.vcf_anno_exome}

		# Annotate with genome data
		bcftools annotate \
			-a {params.gnomad_genome} \
			-h {params.gnomad_genome_header} \
			-c CHROM,POS,POS,REF,ALT,AF_genome,AF_popmax_genome \
			{output.vcf_anno_exome} \
			-o {output.vcf_anno_genome} \
			-Oz >> {log} 2>&1
		tabix -p vcf {output.vcf_anno_genome}
		'''

rule process_tier:
	input:
		vcf_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz.tbi",
		vcf_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz.tbi",
		vcf_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/08_anno_gnomAD_dbsnp/{model}/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz.tbi"
	output:
		vcf_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
		tbi_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.SNV.{ref_version}.vcf.gz.tbi",
		vcf_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.INS.{ref_version}.vcf.gz",
		tbi_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.INS.{ref_version}.vcf.gz.tbi",
		vcf_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
		tbi_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.DEL.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{cov}.SNV.{ref_version}.tier.logs"
	params:
		process_tier_script=config['process_tier_script']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		python {params.process_tier_script} {input.vcf_snv} {output.vcf_snv}
		tabix -p vcf {output.vcf_snv}

		python {params.process_tier_script} {input.vcf_ins} {output.vcf_ins}
		tabix -p vcf {output.vcf_ins}

		python {params.process_tier_script} {input.vcf_del} {output.vcf_del}
		tabix -p vcf {output.vcf_del}
		'''

rule Annotation_SNV:
	input:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.SNV.{ref_version}.vcf.gz.tbi"
	output:
		outputanno="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt"
	params:
		annovar_dir=config['annovar_dir'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		outputanno="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoSNV"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		if [ $(zgrep -v '^#' {input.vcf} | wc -l) != 0 ]; then
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
		else
		 # VCF is empty, create empty output file
			touch {output.outputanno}
			touch {params.outputanno}.{ref_version}_multianno.vcf
		fi
		'''

rule Annotation_INS:
	input:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.INS.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.INS.{ref_version}.vcf.gz.tbi"
	output:
		outputanno="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoINS.{ref_version}_multianno.txt"
	params:
		annovar_dir=config['annovar_dir'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		outputanno="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoINS"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		if [ $(zgrep -v '^#' {input.vcf} | wc -l) != 0 ]; then
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
		else
		 # VCF is empty, create empty output file
			touch {output.outputanno}
			touch {params.outputanno}.{ref_version}_multianno.vcf
		fi
		'''

rule Annotation_DEL:
	input:
		vcf="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/01_mutect2/02_mosaicforecast/09_mosaic_tier/{model}/{sample}.{cov}.DEL.{ref_version}.vcf.gz.tbi"
	output:
		outputanno="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
	params:
		annovar_dir=config['annovar_dir'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		outputanno="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoDEL"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		# Check if VCF has variants (non-header lines)
		# VCF has variants, run ANNOVAR
		if [ $(zgrep -v '^#' {input.vcf} | wc -l) != 0 ]; then
			perl {params.annovar_dir}/table_annovar.pl \
			{input.vcf} \
			{params.annovar_dir}/humandb_{ref_version} \
			-buildver {ref_version} \
			-out {params.outputanno} \
			-remove \
			-protocol refGene,dbnsfp42a,clinvar_20240917,gnomad41_genome,gnomad41_exome \
			-operation g,f,f \
			-nastring . \
			-vcfinput
		else
		 # VCF is empty, create empty output file
			touch {output.outputanno}
			touch {params.outputanno}.{ref_version}_multianno.vcf
		fi
		'''

rule reformat_annotation_clinvar:
	input:
		anno_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		anno_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoINS.{ref_version}_multianno.txt",
		anno_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
	output:
		anno_snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		anno_ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.annoINS.{ref_version}_multianno.txt",
		anno_del="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/04_deepvariant/logs/{sample}.{cov}.{ref_version}.m2.reforma_anno.log"
	params:
		reformat_script=config['reformat_script'],
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		python {params.reformat_script} {input.anno_snv} {output.anno_snv} --info-column "Otherinfo12" --gt-column "Otherinfo13" > {log} 2>&1
		python {params.reformat_script} {input.anno_ins} {output.anno_ins} --info-column "Otherinfo12" --gt-column "Otherinfo13" >> {log} 2>&1
		python {params.reformat_script} {input.anno_del} {output.anno_del} --info-column "Otherinfo12" --gt-column "Otherinfo13" >> {log} 2>&1
		"""

rule annotate_rcnv_gnomadlof_snv:
	input:
		tier_anno="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt"
	output:
		txt="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annocnv_gnomadlof.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{cov}.{ref_version}.m2.rcnv_lof.log"
	params:
		ref_version=config['ref_version'],
		gnomad_LoF=config['gnomad_LoF'],
		rCNV_gene_score=config['rCNV_gene_score'],
		sub="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.SNV.{ref_version}.exonic_splicing_multianno.txt",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		"""
		cat <(awk '{{if($6=="exonic"){{print $0}}}}' {input.tier_anno} | grep -E 'nonsynonymous|stop') <(awk '{{if($6=="splicing"){{print $0}}}}' {input.tier_anno}) > {params.sub}
		cat <(paste <(head -n 1 {input.tier_anno}) <(head -n 1 {params.gnomad_LoF}) <(head -n 1 {params.rCNV_gene_score})) <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$7]){{print $0"\t"c[$7]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.gnomad_LoF} <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$7]){{print $0"\t"c[$7]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.rCNV_gene_score} {params.sub})) > {output.txt}
		rm {params.sub}
		"""

rule reformat_rcnv_gnomadlof:
	input:
		txt="{outpath}/03_variants/01_mutect2/02_mosaicforecast/10_annovar/{model}/{sample}.{cov}.annocnv_gnomadlof.{ref_version}_multianno.txt"
	output:
		txt="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.exonic.splicing.{ref_version}.txt"
	log:
		"{outpath}/03_variants/04_deepvariant/logs/{sample}.{cov}.{ref_version}.m2.reformat.log"
	params:
		reformat_script=config['reformat_script'],
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		python {params.reformat_script} {input.txt} {output.txt} --info-column "Otherinfo12" --gt-column "Otherinfo13" > {log} 2>&1
		"""

rule summary:
	input:
		snv_exonic="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.exonic.splicing.{ref_version}.txt",
		snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.annoINS.{ref_version}_multianno.txt",
		deletion="{outpath}/03_variants/01_mutect2/02_mosaicforecast/11_annovar_reformat/{model}/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
	output:
		snv="{outpath}/03_variants/01_mutect2/02_mosaicforecast/{model}/{sample}.{cov}.SNV.{ref_version}.mosaic_summary.txt",
		ins="{outpath}/03_variants/01_mutect2/02_mosaicforecast/{model}/{sample}.{cov}.INS.{ref_version}.mosaic_summary.txt",
		deletion="{outpath}/03_variants/01_mutect2/02_mosaicforecast/{model}/{sample}.{cov}.DEL.{ref_version}.mosaic_summary.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		wc -l {input.snv} > {output.snv}
		wc -l {input.ins} > {output.ins}
		wc -l {input.deletion} > {output.deletion}
		''' 