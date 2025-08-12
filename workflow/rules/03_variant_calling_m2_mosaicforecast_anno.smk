
rule MT2_initial_filter:
	input:
		vcf="{outpath}/03_variants/mutect2/00_initial_call/03_filter_mutect_call/{sample}.vcf.gz"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/02_init_filter/{sample}.mt2pon.AF0.02.bed"
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
		mt2pon="{outpath}/03_variants/mutect2/00_mosaicforecast/02_init_filter/{sample}.mt2pon.AF0.02.bed"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.noSegDup.bed"
	params:
		segdup=config['SegDup_and_clustered']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	conda:
		config["bcftools_1.9"]
	shell:
		'''
		subtractBed -a {input.mt2pon} -b {params.segdup} > {output}
		'''

rule annovar_formatter:
	input:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.noSegDup.bed"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
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
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.SNV.chr.bed"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)==1 && length($5)==1'\'' |awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"


rule MAF0_extraction_INS:
	input:
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.INS.chr.bed"
	params:
		allrepeats_forindel=config['allrepeats_forindel']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	conda:
		config["bcftools_1.9"]
	shell:
		"subtractBed -a <(cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' |awk '\''length($4)< length($5)'\'') -b {params.allrepeats_forindel} | awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"

rule MAF0_extraction_DEL:
	input:
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/03_repeat_filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.DEL.chr.bed"
	params:
		allrepeats_forindel=config['allrepeats_forindel']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	conda:
		config["bcftools_1.9"]
	shell:
		"subtractBed -a <(cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)> length($5)'\'') -b {params.allrepeats_forindel} | awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"

rule split_bed:
	input:
		bed_snv="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.SNV.chr.bed",
		bed_ins="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.INS.chr.bed",
		bed_del="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.DEL.chr.bed"
	output:
		bed_chunks_snv="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.snv.bed",
		bed_chunks_ins="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.ins.bed",
		bed_chunks_del="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.del.bed"
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
		bed="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL_split/{sample}.{individual_chr}.snv.bed",
		bam="{outpath}/02_map/08_bqsr/{sample}.bam",
		bai="{outpath}/02_map/08_bqsr/{sample}.bai"
	output:
		feature="{outpath}/03_variants/mutect2/00_mosaicforecast/05_feature_each_chromosome/{sample}.SNV.{individual_chr}.features"
	params:
		outdir="{outpath}/02_map/08_bqsr",
		umap=config['umap'],
		ReadLevel_Features_extraction=config['ReadLevel_Features_extraction'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/{sample}.SNV.{individual_chr}.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["terra_py_tools"]
	shell:
		'''
		export PYTHONWARNINGS="ignore" 
		python {params.ReadLevel_Features_extraction} \
		{input.bed} \
		{output.feature} \
		{params.outdir} \
		{params.ref} \
		{params.umap} \
		1 \
		bam > {log} 2>&1 
		'''

rule merge_features_SNV:
	input:
		lambda wildcards: [f"{wildcards.outpath}/03_variants/mutect2/00_mosaicforecast/05_feature_each_chromosome/{wildcards.sample}.SNV.{chr}.features" for chr in CHROMOSOMES]
		#expand("{{outpath}}/03_variants/mutect2/00_mosaicforecast/05_feature_each_chromosome/{{sample}}.SNV.{chr}.features" for chr in CHROMOSOMES)
	output:
		feature="{outpath}/03_variants/mutect2/00_mosaicforecast/05_feature/{sample}.SNV.features"
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
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.INS.chr.bed",
		file2="{outpath}/02_map/08_bqsr/{sample}.bam"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/05_feature/{sample}.INS.features"
	params:
		outdir="{outpath}/02_map/08_bqsr/",
		umap=config['umap'],
		ReadLevel_Features_extraction=config['ReadLevel_Features_extraction'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/{sample}.INS.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["terra_py_tools"]
	shell:
		'''
		export PYTHONWARNINGS="ignore" 
		python {params.ReadLevel_Features_extraction} \
		{input.file1} \
		{output} \
		{params.outdir} \
		{params.ref} \
		{params.umap} 1 bam > {log} 2>&1
		'''

rule feature_extraction_DEL:
	input:
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/04_extrac_SNV_INDEL/{sample}.MAF0.DEL.chr.bed",
		file2="{outpath}/02_map/08_bqsr/{sample}.bam"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/05_feature/{sample}.DEL.features"
	params:
		outdir="{outpath}/02_map/08_bqsr/",
		umap=config['umap'],
		ReadLevel_Features_extraction=config['ReadLevel_Features_extraction'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/{sample}.DEL.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["terra_py_tools"]
	shell:
		'''
		export PYTHONWARNINGS="ignore" 
		python {params.ReadLevel_Features_extraction} \
		{input.file1} \
		{output} \
		{params.outdir} \
		{params.ref} \
		{params.umap} \
		1 \
		bam > {log} 2>&1
		'''

rule g1000_avail_acess_filter_SNV:
	input:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/05_feature/{sample}.SNV.features"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/06_1000g_avail_access/{sample}.SNV.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		'''
		cat <(head -n 1 {input}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input} | cut -f 2-3) {input}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_INS:
	input:
		feature="{outpath}/03_variants/mutect2/00_mosaicforecast/05_feature/{sample}.INS.features"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/06_1000g_avail_access/{sample}.INS.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		'''
		cat <(head -n 1 {input.feature}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_DEL:
	input:
		feature="{outpath}/03_variants/mutect2/00_mosaicforecast/05_feature/{sample}.DEL.features"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/06_1000g_avail_access/{sample}.DEL.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["bcftools_1.9"]
	shell:
		'''
		cat <(head -n 1 {input.feature}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		'''

rule Prediction_SNV:
	input:
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/06_1000g_avail_access/{sample}.SNV.no1000g.features"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.SNV.predictions"
	params:
		refine_beta=config['refine_beta'] if config['mosaic_pred_model'] == 'Refine' else config['phase_model'],
		model=config['mosaic_pred_model'],
		prediction_r_script=config['prediction_r_script']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["terra_r"]
	shell:
		'''
		Rscript {params.prediction_r_script} {input.file1} {params.refine_beta} {params.model} {output}
		'''

rule Prediction_INS:
	input:
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/06_1000g_avail_access/{sample}.INS.no1000g.features"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.INS.predictions"
	params:
		model=config['mosaic_pred_model'],
		refine_beta=config['insertion_model'],
		prediction_r_script=config['prediction_r_script']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["terra_r"]
	shell:
		'''
		Rscript {params.prediction_r_script} {input.file1} {params.refine_beta} {params.model} {output}
		'''

rule Prediction_DEL:
	input:
		file1="{outpath}/03_variants/mutect2/00_mosaicforecast/06_1000g_avail_access/{sample}.DEL.no1000g.features"
	output:
		"{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.DEL.predictions"
	params:
		model=config['mosaic_pred_model'],
		refine_beta=config['deletion_model'],
		prediction_r_script=config['prediction_r_script']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["terra_r"]
	shell:
		'''
		Rscript {params.prediction_r_script} {input.file1} {params.refine_beta} {params.model} {output}
		'''

rule extract_bed_snv:
	input:
		prediction_snv="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.SNV.predictions"
	output:
		bed_snv="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.SNV.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.SNV.{ref_version}.bed.temp_file",
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
		prediction_ins="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.INS.predictions"
	output:
		bed_ins="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.INS.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.INS.{ref_version}.bed.temp_file"
	shell:
		'''
		awk '$35=="hap=3"{{print $0}}' {input.prediction_ins} > {params.temp_file}
		if [ -s {params.temp_file} ]; then
			cat {params.temp_file} | cut -f 1 | cut -d"~" -f2-5 | sed 's/~/\t/g' | grep -v "id" | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n > {output.bed_ins}
		else
			: > {output.bed_ins}
		fi
		rm {params.temp_file}
		'''

rule extract_bed_del:
	input:
		prediction_del="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.DEL.predictions"
	output:
		bed_del="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.DEL.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.DEL.{ref_version}.bed.temp_file"
	shell:
		'''
		awk '$35=="hap=3"{{print $0}}' {input.prediction_del} > {params.temp_file}
		if [ -s {params.temp_file} ]; then
			cat {params.temp_file} | cut -f 1 | cut -d"~" -f2-5 | sed 's/~/\t/g' | grep -v "id" | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n > {output.bed_del}
		else
			: > {output.bed_del}
		fi
		rm {params.temp_file}
		'''

rule merge_bed_snv_indel:
	input:
		bed_del="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.DEL.{ref_version}.bed",
		bed_ins="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.INS.{ref_version}.bed",
		bed_snv="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.SNV.{ref_version}.bed"
	output:
		bed_all="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.all.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']	
	shell:
		'''
		cat {input.bed_snv} {input.bed_del} {input.bed_ins} > {output.bed_all}
		'''


rule extrac_all_subvcf:
	input:
		vcf="{outpath}/03_variants/mutect2/00_initial_call/03_filter_mutect_call/{sample}.mt2pon.filter.vcf.gz",
		bed="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.all.{ref_version}.bed"
	output:
		vcf="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.all.{ref_version}.vcf.gz",
		dpaf="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.all.{ref_version}.Geno.DP.AF.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config["bedtools_2.31.1"]
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

rule anno_gnomAD_dbsbp:
	input:
		vcf="{outpath}/03_variants/mutect2/00_mosaicforecast/07_mosaic_prediction/{sample}.{model}.{cov}.all.{ref_version}.vcf.gz",
	output:
		vcf_anno_dbsnp="{outpath}/03_variants/mutect2/00_mosaicforecast/08_anno_gnomAD_dbsnp/{sample}.{model}.{cov}.all.anno.dbsnp.{ref_version}.vcf.gz",
		vcf_anno_exome="{outpath}/03_variants/mutect2/00_mosaicforecast/08_anno_gnomAD_dbsnp/{sample}.{model}.{cov}.all.anno.exome.{ref_version}.vcf.gz",
		vcf_anno_genome="{outpath}/03_variants/mutect2/00_mosaicforecast/08_anno_gnomAD_dbsnp/{sample}.{model}.{cov}.all.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_anno_genome="{outpath}/03_variants/mutect2/00_mosaicforecast/08_anno_gnomAD_dbsnp/{sample}.{model}.{cov}.all.anno.exome.genome.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{model}.{cov}.all.{ref_version}.dbsnp.gnomad.logs"
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
	container:
		config["bcftools_1.9"]
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
		vcf="{outpath}/03_variants/mutect2/00_mosaicforecast/08_anno_gnomAD_dbsnp/{sample}.{model}.{cov}.all.anno.exome.genome.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/mutect2/00_mosaicforecast/08_anno_gnomAD_dbsnp/{sample}.{model}.{cov}.all.anno.exome.genome.{ref_version}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/mutect2/00_mosaicforecast/09_mosaic_tier/{sample}.{model}.{cov}.all.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/mutect2/00_mosaicforecast/09_mosaic_tier/{sample}.{model}.{cov}.all.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{model}.{cov}.all.{ref_version}.tier.logs"
	params:
		process_tier_script=config['process_tier_script']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["python_alpine3.21"]
	shell:
		'''
		python {params.process_tier_script} {input.vcf} {output.vcf}
		tabix -p vcf {output.vcf}
		'''

rule Annotation:
	input:
		vcf="{outpath}/03_variants/mutect2/00_mosaicforecast/09_mosaic_tier/{sample}.{model}.{cov}.all.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/mutect2/00_mosaicforecast/09_mosaic_tier/{sample}.{model}.{cov}.all.{ref_version}.vcf.gz.tbi"
	output:
		outputanno="{outpath}/03_variants/mutect2/00_mosaicforecast/10_annovar/{sample}.{model}.{cov}.all.{ref_version}_multianno.txt"
	params:
		annovar_dir=config['annovar_dir'],
		ref_version=config['ref_version'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		outputanno="{outpath}/03_variants/mutect2/00_mosaicforecast/10_annovar/{sample}.{model}.{cov}.all"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["terra_perl_anno"]
	shell:
		'''
		if [ $(zgrep -v '^#' {input.vcf} | wc -l) != 0 ]; then
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
		else
		 # VCF is empty, create empty output file
			touch {output.outputanno}
			touch {params.outputanno}.{params.ref_version}_multianno.vcf
		fi
		'''

rule reformat_annotation_clinvar:
	input:
		anno="{outpath}/03_variants/mutect2/00_mosaicforecast/10_annovar/{sample}.{model}.{cov}.all.{ref_version}_multianno.txt"
	output:
		anno="{outpath}/03_variants/mutect2/00_mosaicforecast/11_annovar_reformat/{sample}.{model}.{cov}.all.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/04_deepvariant/logs/{sample}.{model}.{cov}.{ref_version}.m2.reforma_anno.log"
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
		python {params.reformat_script} {input.anno} {output.anno} --info-column "Otherinfo12" --gt-column "Otherinfo13" > {log} 2>&1
		"""

rule annotate_rcnv_gnomadlof:
	input:
		tier_anno="{outpath}/03_variants/mutect2/00_mosaicforecast/10_annovar/{sample}.{model}.{cov}.all.{ref_version}_multianno.txt"
	output:
		txt="{outpath}/03_variants/mutect2/00_mosaicforecast/10_annovar/{sample}.{model}.{cov}.annocnv_gnomadlof.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{model}.{cov}.{ref_version}.m2.rcnv_lof.log"
	params:
		ref_version=config['ref_version'],
		gnomad_LoF=config['gnomad_LoF'],
		rCNV_gene_score=config['rCNV_gene_score'],
		sub="{outpath}/03_variants/mutect2/00_mosaicforecast/10_annovar/{sample}.{model}.{cov}.SNV.{ref_version}.exonic_splicing_multianno.txt",
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
		txt="{outpath}/03_variants/mutect2/00_mosaicforecast/10_annovar/{sample}.{model}.{cov}.annocnv_gnomadlof.{ref_version}_multianno.txt"
	output:
		txt="{outpath}/03_variants/mutect2/00_mosaicforecast/11_annovar_reformat/{sample}.{model}.{cov}.exonic.splicing.{ref_version}.txt"
	log:
		"{outpath}/03_variants/04_deepvariant/logs/{sample}.{model}.{cov}.{ref_version}.m2.reformat.log"
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

rule summary:
	input:
		ins=expand([
			"{outpath}/03_variants/mutect2/00_mosaicforecast/11_annovar_reformat/{sample}.{model}.{cov}.all.{ref_version}_multianno.txt",
			"{outpath}/03_variants/mutect2/00_mosaicforecast/11_annovar_reformat/{sample}.{model}.{cov}.exonic.splicing.{ref_version}.txt"
			],
			outpath=Outpath,
			model=Model,
			cov=Coverage,
			sample=Sample,
			ref_version=Ref_version)
	output:
		snv="{outpath}/03_variants/mosaic_report.html"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		echo {output.snv}
		'''