rule Mutect2:
	message:
		"""
		---
		using mutect2 4.6.1.0
		---
		"""
	input:
		"{outpath}/02_map/bqsr/{sample}/{sample}.bam"
	output:
		raw_vcf="{outpath}/03_variants/mutect2/result/sub/{sample}/{sample}.{individual_chr}.mt2pon.vcf.gz",
		raw_stat="{outpath}/03_variants/mutect2/result/sub/{sample}/{sample}.{individual_chr}.mt2pon.vcf.gz.stats"
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
		lambda wildcards: [f"{wildcards.outpath}/03_variants/mutect2/result/sub/{wildcards.sample}/{wildcards.sample}.{chr}.mt2pon.vcf.gz" for chr in chromosomes]
	output:
		args="{outpath}/03_variants/mutect2/result/{sample}/{sample}.args",
		vcf="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.vcf.gz",
		idx="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.vcf.gz.tbi"
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
		lambda wildcards: [f"{wildcards.outpath}/03_variants/mutect2/result/sub/{wildcards.sample}/{wildcards.sample}.{chr}.mt2pon.vcf.gz.stats" for chr in chromosomes]
		#expand("{outpath}/{key}/03_variants/mutect2/result/sub/{key}.{chr}.mt2pon.vcf.gz.stats", outpath=outpath, key=sample.keys(), chr=chromosomes)
	output:
		"{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.vcf.gz.stats"
	params:
		list_para = lambda wildcards: ' '.join([f"--stats {wildcards.outpath}/03_variants/mutect2/result/sub/{wildcards.sample}/{wildcards.sample}.{chr}.mt2pon.vcf.gz.stats" for chr in chromosomes]),
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
		java -Xms1g -XX:ParallelGCThreads={threads} -jar {params.gatk} MergeMutectStats {params.list_para} -O {output} > {log} 2>&1
		'''

rule FilterMutectCall:
	input:
		vcf="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.vcf.gz",
		stats="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.vcf.gz.stats",
		index="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.vcf.gz.tbi"
	output:
		"{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.filter.vcf.gz"
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
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} FilterMutectCalls -R {params.ref} -V {input.vcf} -O {output} > {log} 2>&1
		'''

rule MT2_initial_filter:
	input:
		vcf="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.filter.vcf.gz"
	output:
		"{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.bed"
	params:
		sample="{sample}",
		af_threshold="0.03" if config["pcr_based"] == True else "0.02"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat <(zcat {input.vcf} | grep -v "^#"| grep PASS | gawk '{{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){{print $0}}}}' | cut -f1,2,4,5,10 | sed "s/:/\\t/g" | sed "s/,/\\t/g" | awk '$8>=0.03 && $8<0.4' | grep -v '0|1' | grep -v '1|0') <(zcat {input.vcf} | grep -v "^#" | grep PASS | gawk '{{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){{print $0}}}}' | cut -f1,2,4,5,10 | sed "s/:/\\t/g" | sed "s/,/\\t/g" | awk '$8>={params.af_threshold} && $8<0.4 ' | grep -E "0\|1|1\|0") | cut -f 1-4,6-8 | awk '{{OFS="\\t";print $1,$2-1,$2,$3,$4,\"{params.sample}\",$5,$6,$7}}' > {output}
		'''

rule repeat_filter:
	input:
		mt2pon="{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.bed"
	output:
		"{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.noSegDup.bed"
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
		"{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.noSegDup.bed"
	output:
		"{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.ANNOVAR.list"
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
		file1="{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/mutect2/filter/{sample}/{sample}.MAF0.SNV.chr.bed"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)==1 && length($5)==1'\'' |awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"


rule MAF0_extraction_INS:
	input:
		file1="{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/mutect2/filter/{sample}/{sample}.MAF0.INS.chr.bed"
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
		file1="{outpath}/03_variants/mutect2/filter/{sample}/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/mutect2/filter/{sample}/{sample}.MAF0.DEL.chr.bed"
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

rule feature_extraction_SNV:
	input:
		file1="{outpath}/03_variants/mutect2/filter/{sample}/{sample}.MAF0.SNV.chr.bed",
		file2="{outpath}/02_map/bqsr/{sample}/{sample}.bam"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.SNV.features"
	params:
		outdir="{outpath}/02_map/bqsr/{sample}/",
		umap=config['umap'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/{sample}.SNV.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate py3.7.1
		export PYTHONWARNINGS="ignore" 
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py {input.file1} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1 
		conda deactivate
		'''

rule feature_extraction_INS:
	message:
		"since snakemake cannot create a conda env from ../envs/py3.7.1.yaml, but we would liek to use --use-conda, so we use conda activate"
	input:
		file1="{outpath}/03_variants/mutect2/filter/{sample}/{sample}.MAF0.INS.chr.bed",
		file2="{outpath}/02_map/bqsr/{sample}/{sample}.bam"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.INS.features"
	params:
		outdir="{outpath}/02_map/bqsr/{sample}/",
		umap=config['umap'],
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
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py {input.file1} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1
		conda deactivate
		'''

rule feature_extraction_DEL:
	input:
		file1="{outpath}/03_variants/mutect2/filter/{sample}/{sample}.MAF0.DEL.chr.bed",
		file2="{outpath}/02_map/bqsr/{sample}/{sample}.bam"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.DEL.features"
	params:
		outdir="{outpath}/02_map/bqsr/{sample}/",
		umap=config['umap'],
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
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py \
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
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.SNV.features"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.SNV.no1000g.features"
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
		feature="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.INS.features"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.INS.no1000g.features"
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
		feature="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.DEL.features"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.DEL.no1000g.features"
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
		file1="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.SNV.no1000g.features"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.predictions"
	params:
		refine_beta=config['refine_beta']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Refine {output}
		'''

rule Prediction_INS:
	input:
		file1="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.INS.no1000g.features"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Refine {output}
		'''

rule Prediction_DEL:
	input:
		file1="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.DEL.no1000g.features"
	output:
		"{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Refine {output}
		'''

rule Annotation_SNV_Mosaic:
	input:
		file1="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.predictions"
	output:
		inputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.input.predictions.txt",
		outputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/03_variants/mutect2/feature/{sample}",
		outputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.annoSNV"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		grep "mosaic" {input.file1} | cut -f 1,35 | awk 'BEGIN {{ FS = "~" }} ; {{ print $2"\\t"$3"\\t"$3"\\t"$4"\\t"$5"\\t"$6 }}' | tail -n +2 | awk '$6=="mosaic" {{print $0}}' > {output.inputanno}
		cd {params.outdir}
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl {output.inputanno} /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} -buildver {ref_version} -out {params.outputanno} -remove -protocol refGene,dbnsfp42a -operation g,f -nastring . 
		'''

rule Annotation_INS_Mosaic:
	input:
		file1="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.predictions"
	output:
		inputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.{ref_version}.input.predictions.txt",
		mosaic="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.mosaic.predictions.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/03_variants/mutect2/feature/{sample}",
		outputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.ouput.predictions"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		#Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/mosaic_bin/choooseMosaicVar.R -i {input.file1} -o {output.mosaic} -t {output.inputanno}
		grep "hap=3" {input.file1} > {output.mosaic}
		cat {output.mosaic} | awk 'BEGIN {{ FS = "~" }} ; {{ print $2"\\t"$3"\\t"$3"\\t"$4"\\t"$5 }}' | tail -n +2 > {output.inputanno}
		cd {params.outdir}
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl {output.inputanno} /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} -buildver {ref_version} -out {params.outputanno} -remove -protocol refGene,dbnsfp42a -operation g,f -nastring . 
		'''

rule Annotation_DEL_Mosaic:
	input:
		file1="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.predictions"
	output:
		inputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.{ref_version}.input.predictions.txt",
		mosaic="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.mosaic.predictions.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/03_variants/mutect2/feature/{sample}",
		outputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.ouput.predictions"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		grep "hap=3" {input.file1} > {output.mosaic}
		cat {output.mosaic} | awk 'BEGIN {{ FS = "~" }} ; {{ print $2"\\t"$3"\\t"$3+length($4)-1"\\t"$4"\\t"$5 }}' | tail -n +2 > {output.inputanno}
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl {output.inputanno} /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} -buildver {ref_version} -out {params.outputanno} -remove -protocol refGene,dbnsfp42a -operation g,f -nastring . 
		'''

rule extrac_SNVsubvcf:
	input:
		vcf="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.filter.vcf.gz",
		inputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf",
		tmp="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf.tmp",
		bed="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.subvcf.bed",
		dpaf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.Geno.DP.AF.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		awk '$6=="mosaic"{{print $0}}' {input.inputanno} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2}}' | sort -k1,1 -k2,2n > {output.bed}
		tabix {input.vcf} -R {output.bed} > {output.tmp}
		cat <(zcat {input.vcf} | grep "^#") <(cat {output.tmp}) > {output.vcf}
		cut -f 1,2,4,5,10 {output.tmp} | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		'''

rule extrac_INSsubvcf:
	input:
		vcf="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.filter.vcf.gz",
		inputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.{ref_version}.sub.vcf",
		tmp="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.{ref_version}.sub.vcf.tmp",
		bed="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.{ref_version}.subvcf.bed",
		dpaf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.{ref_version}.Geno.DP.AF.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat {input.inputanno} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2}}' | sort -k1,1 -k2,2n > {output.bed}
		tabix {input.vcf} -R {output.bed} > {output.tmp}
		cat <(zcat {input.vcf} | grep "^#") <(cat {output.tmp}) > {output.vcf}
		cut -f 1,2,4,5,10 {output.tmp} | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		'''
		
rule extrac_DELsubvcf:
	input:
		vcf="{outpath}/03_variants/mutect2/result/{sample}/{sample}.mt2pon.filter.vcf.gz",
		inputanno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.{ref_version}.sub.vcf",
		tmp="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.{ref_version}.sub.vcf.tmp",
		bed="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.{ref_version}.subvcf.bed",
		dpaf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.{ref_version}.Geno.DP.AF.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat {input.inputanno} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2}}' | sort -k1,1 -k2,2n > {output.bed}
		tabix {input.vcf} -R {output.bed} > {output.tmp}
		cat <(zcat {input.vcf} | grep "^#") <(cat {output.tmp}) > {output.vcf}
		cut -f 1,2,4,5,10 {output.tmp} | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		'''

rule gnomAD_vcf_tier2_3_subchr:
	input:
		vcf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf",
	output:
		vcf_tier3="{outpath}/03_variants/mutect2/tier/{sample}/overlap_gnomAD/{sample}.{cov}.SNV.{ref_version}.sub.tier3.{individual_chr}.vcf",
		vcf_tier2="{outpath}/03_variants/mutect2/tier/{sample}/overlap_gnomAD/{sample}.{cov}.SNV.{ref_version}.sub.tier2.{individual_chr}.vcf"
	params:
		gnomad_2_1_1_0_001="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_gnomad211_genome_exome/sub_chr/{ref_version}_gnomad211_exome_genome_afpop_gt_0.001_chr_snponly_{individual_chr}.txt",
		gnomad_2_1_1_0="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_gnomad211_genome_exome/sub_chr/{ref_version}_gnomad211_exome_genome_afpop_gt_0_le_0.001_chr_snponly_{individual_chr}.txt"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		awk 'NR==FNR{{c[$1,$2,$4,$5]=$0}}NR!=FNR{{if(c[$1,$2,$4,$5]) {{print $0}}}}' {params.gnomad_2_1_1_0_001} {input.vcf} > {output.vcf_tier3}
		awk 'NR==FNR{{c[$1,$2,$4,$5]=$0}}NR!=FNR{{if(c[$1,$2,$4,$5]) {{print $0}}}}' {params.gnomad_2_1_1_0} {input.vcf} > {output.vcf_tier2}
		'''

rule gather_gnomAD_vcf_tier2_3:
	input:
		vcf_tier3=lambda wildcards: [f"{wildcards.outpath}/03_variants/mutect2/tier/{wildcards.sample}/overlap_gnomAD/{wildcards.sample}.{wildcards.cov}.SNV.{wildcards.ref_version}.sub.tier3.{chr}.vcf" for chr in chromosomes],
		vcf_tier2=lambda wildcards: [f"{wildcards.outpath}/03_variants/mutect2/tier/{wildcards.sample}/overlap_gnomAD/{wildcards.sample}.{wildcards.cov}.SNV.{wildcards.ref_version}.sub.tier2.{chr}.vcf" for chr in chromosomes]		
	output:
		vcf_tier3="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier3.vcf",
		vcf_tier2="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier2.vcf"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		cat {input.vcf_tier3} > {output.vcf_tier3}
		cat {input.vcf_tier2} > {output.vcf_tier2}
		'''

rule gnomAD_vcf_tier1:
	input:
		vcf="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf",
		vcf_tier3="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier3.vcf",
		vcf_tier2="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier2.vcf"
	output:
		vcf_tier1="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier1.vcf"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		awk 'NR==FNR{{c[$1,$2,$4,$5]=$0}}NR!=FNR{{if(!c[$1,$2,$4,$5]) {{print $0}}}}' <(cat {input.vcf_tier2} {input.vcf_tier3}) {input.vcf} > {output.vcf_tier1}
		'''

rule mark_tier123:
	input:
		anno="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		vcf_tier3="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier3.vcf",
		vcf_tier2="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier2.vcf",
		vcf_tier1="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier1.vcf"
	output:
		tier_anno="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.tier.{ref_version}_multianno.txt"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		cat <(awk 'BEGIN{{OFS="\t"}}NR==FNR{{c[$1,$2,$4,$5]=$0}} NR!=FNR{{if (c[$1,$2,$4,$5]) {{print "tier1\t"$0}}}}' {input.vcf_tier1} {input.anno}) <(awk 'BEGIN{{OFS="\t"}}NR==FNR{{c[$1,$2,$4,$5]=$0}} NR!=FNR{{if (c[$1,$2,$4,$5]) {{print "tier2\t"$0}}}}' {input.vcf_tier2} {input.anno}) <(awk 'BEGIN{{OFS="\t"}}NR==FNR{{c[$1,$2,$4,$5]=$0}} NR!=FNR{{if (c[$1,$2,$4,$5]) {{print "tier3\t"$0}}}}' {input.vcf_tier3} {input.anno}) > {output.tier_anno}
		'''

rule summary:
	input:
		#if another job is running, same output "{outpath}/mosaic_summary.txt" will be rewriten or will stop in advance.
		#expand(["{outpath}/03_variants/mutect2/tier/{u.sample}/{u.sample}.{cov}.SNV.tier.{ref_version}_multianno.txt"], outpath=outpath, u=units.itertuples(), cov=cov, ref_version=ref_version)
		snv="{outpath}/03_variants/mutect2/tier/{sample}/{sample}.{cov}.SNV.tier.{ref_version}_multianno.txt",
		ins="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.INS.{ref_version}.Geno.DP.AF.txt",
		del1="{outpath}/03_variants/mutect2/feature/{sample}/{sample}.{cov}.DEL.{ref_version}.Geno.DP.AF.txt"
	output:
		snv="{outpath}/{sample}.{cov}.SNV.tier.{ref_version}.mosaic_summary.txt",
		ins="{outpath}/{sample}.{cov}.INS.{ref_version}.mosaic_summary.txt",
		del1="{outpath}/{sample}.{cov}.DEL.{ref_version}.mosaic_summary.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		wc -l {input.snv} > {output.snv}
		wc -l {input.ins} > {output.ins}
		wc -l {input.del1} > {output.del1}
		'''
