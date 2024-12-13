rule Mutect2:
	message:
		"""
		---
		using mutect2 4.0.12.0 instead of 4.1.8.1 to call variants, since 4.1.8.1 didnt take account of overlapping of paired end reads.
		---
		"""
	input:
		"{outpath}/{sample}/02_map/bqsr/{sample}.sort.rmdup.bqsr.bam"
	output:
		raw_vcf="{outpath}/{sample}/03_variants/mutect2/result/sub/{sample}.{individual_chr}.MT2PON.vcf.gz",
		raw_stat="{outpath}/{sample}/03_variants/mutect2/result/sub/{sample}.{individual_chr}.MT2PON.vcf.gz.stat"
	log:
		"{outpath}/{sample}/03_variants/mutect2/logs/{sample}.{individual_chr}.mutect.log"
	threads:
		8
	params:
		mem="4000",
		pon=config['pon'],
		af_only_gnomad=config['af_only_gnomad'],
		ref=config['reference'],
		sample="{sample}",
		interval_list=intervals_dir + "/{individual_chr}.intervals.list",
		gatk=config['gatk_current_using']
	shell:
		'''
		set +u
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.0.xx 
		java -Xms10g -XX:ParallelGCThreads={threads} -jar {params.gatk} Mutect2 -R {params.ref} -I {input} --pon {params.pon} -tumor {params.sample} --germline-resource {params.af_only_gnomad} -L {params.interval_list} --interval-padding 100 -O {output.raw_vcf} > {log} 2>&1
		conda deactivate
		set -u
		'''

rule vcf_args:
	message:
		"""
		---
		merge the vcf file for individual chromosome into single one vcf for each sample
		---
		"""
	input:
		lambda wildcards: [f"{outpath}/{wildcards.sample}/03_variants/mutect2/result/sub/{wildcards.sample}.{chr}.MT2PON.vcf.gz" for chr in chromosomes]
	output:
		"{outpath}/{sample}/03_variants/mutect2/result/{sample}.args"
	log:
		"{outpath}/{sample}/03_variants/mutect2/logs/{sample}.args.logs"
	shell:
		'''
		ls {input} > {output}
		'''

rule vcf_concat:
	message:
		"""
		--- 
		Concat all the sub file to one vcf
		---
		input: {input}
		output: {output}
		"""
	input:
		"{outpath}/{sample}/03_variants/mutect2/result/{sample}.args"
	output:
		"{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.vcf.gz"
	log:
		"{outpath}/{u.sample}/03_variants/mutect2/logs/{sample}.concat.log"
	shell:
		'''
		vcf-concat -f {input}  |  bgzip > {output}
		'''

rule vcf_index:
	message:
		"""
		--- index vcf
		output: {output}
		"""
	input:
		"{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.vcf.gz"
	output:
		protected("{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.vcf.gz.tbi")
	shell:
		'''
		tabix -p vcf {input}
		'''

rule merge_mutectstats:
	input:
		lambda wildcards: [f"{wildcards.outpath}/{wildcards.sample}/03_variants/mutect2/result/sub/{wildcards.sample}.{chr}.MT2PON.vcf.gz.stat" for chr in chromosomes]
		#expand("{outpath}/{key}/03_variants/mutect2/result/sub/{key}.{chr}.MT2PON.vcf.gz.stat", outpath=outpath, key=sample.keys(), chr=chromosomes)
	output:
		"{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.vcf.gz.stats"
	params:
		list_para = lambda wildcards: ' '.join([f"--stats {wildcards.outpath}/{wildcards.sample}/03_variants/mutect2/result/sub/{wildcards.sample}.{chr}.MT2PON.vcf.gz.stat" for chr in chromosomes]),
		gatk=config['gatk_current_using']
	log:
		"{outpath}/{sample}/03_variants/mutect2/logs/{sample}.MergeMutectStats.log"
	shell:
		'''
		java -Xms10g -XX:ParallelGCThreads={threads} -jar {params.gatk} MergeMutectStats {params.list_para} -O {output} > {log} 2>&1
		'''

rule FilterMutectCall:
	input:
		vcf="{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.vcf.gz",
		stats="{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.vcf.gz.stats",
		index="{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.vcf.gz.tbi"
	output:
		"{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.filter.vcf.gz"
	log:
		"{outpath}/{sample}/03_variants/mutect2/logs/{sample}.filter.log"
	threads:
		8
	params:
		mem="4000",
		ref=config['reference'],
		gatk=config['gatk_current_using']
	shell:
		'''
		java -Xms10g -XX:ParallelGCThreads={threads} -jar {params.gatk} FilterMutectCalls -R {params.ref} -V {input.vcf} -O {output} > {log} 2>&1
		'''

rule MT2_initial_filter:
	input:
		vcf="{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.filter.vcf.gz"
	output:
		"{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.bed"
	params:
		sample="{sample}",
		af_threshold="0.03" if config["pcr_based"] == True else "0.02"
	shell:
		'''
		cat <(zcat {input.vcf} | grep -v "^#"| grep PASS | gawk '{{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){{print $0}}}}' | cut -f1,2,4,5,10 | sed "s/:/\\t/g" | sed "s/,/\\t/g" | awk '$8>=0.03 && $8<0.4' | grep -v '0|1' | grep -v '1|0') <(zcat {input.vcf} | grep -v "^#" | grep PASS | gawk '{{match($0,/;POPAF=(([0-9]+\.[0-9]+));/,arr); if(arr[1]!~/-/ && arr[1]>=4){{print $0}}}}' | cut -f1,2,4,5,10 | sed "s/:/\\t/g" | sed "s/,/\\t/g" | awk '$8>={params.af_threshold} && $8<0.4 ' | grep -E "0\|1|1\|0") | cut -f 1-4,6-8 | awk '{{OFS="\\t";print $1,$2-1,$2,$3,$4,\"{params.sample}\",$5,$6,$7}}' > {output}
		'''

rule repeat_filter:
	input:
		mt2pon="{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.bed"
	output:
		"{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.noSegDup.bed"
	params:
		segdup=config['SegDup_and_clustered']
	shell:
		'''
		subtractBed -a {input.mt2pon} -b {params.segdup} > {output}
		'''

rule annovar_formatter:
	input:
		"{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.noSegDup.bed"
	output:
		"{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	shell:
		'''
		cat {input} | awk '{{OFS="\\t\";len=length($4)-length($5);if(len<=0){{print $1,$3,$3,$4,$5,$6}}if(len>0){{print $1,$3,$3+len,$4,$5,$6}}}}'> {output}
		sed -i 's/chr//g' {output}
		'''

rule MAF0_extraction_SNV:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/{sample}/03_variants/mutect2/filter/{sample}.MAF0.SNV.chr.bed"
	shell:
		"cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)==1 && length($5)==1'\'' |awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"


rule MAF0_extraction_INS:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/{sample}/03_variants/mutect2/filter/{sample}.MAF0.INS.chr.bed"
	params:
		allrepeats_forindel=config['allrepeats_forindel']
	shell:
		"subtractBed -a <(cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' |awk '\''length($4)< length($5)'\'') -b {params.allrepeats_forindel} | awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"

rule MAF0_extraction_DEL:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/filter/{sample}.mt2pon.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/{sample}/03_variants/mutect2/filter/{sample}.MAF0.DEL.chr.bed"
	params:
		allrepeats_forindel=config['allrepeats_forindel']
	shell:
		"subtractBed -a <(cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)> length($5)'\'') -b {params.allrepeats_forindel} | awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"

rule rename_bam:
	input:
		bam="{outpath}/{sample}/02_map/bqsr/{sample}.sort.rmdup.bqsr.bam",
		bai="{outpath}/{sample}/02_map/bqsr/{sample}.sort.rmdup.bqsr.bai"
	output:
		bam="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/{sample}.bam",
		bai="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/{sample}.bai"
	params:
		lndir="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/"
	shell:
		'''
			mkdir -p {params.lndir}
			ln -s {input.bam} {output.bam}
			ln -s {input.bai} {output.bai}
		'''

rule feature_extraction_SNV:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/filter/{sample}.MAF0.SNV.chr.bed",
		file2="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/{sample}.bam"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.SNV.features"
	params:
		outdir="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/",
		umap=config['umap'],
		ref=config['reference']
	log:
		"{outpath}/{sample}/03_variants/mutect2/log/{sample}.SNV.features.log"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh && conda activate py3.7.1 
		export PYTHONWARNINGS="ignore" 
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py {input.file1} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1 
		conda deactivate
		'''

rule feature_extraction_INS:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/filter/{sample}.MAF0.INS.chr.bed",
		file2="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/{sample}.bam"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.INS.features"
	params:
		outdir="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/",
		umap=config['umap'],
		ref=config['reference']
	log:
		"{outpath}/{sample}/03_variants/mutect2/log/{sample}.INS.features.log"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh && conda activate py3.7.1 
		export PYTHONWARNINGS="ignore" 
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py {input.file1} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1
		conda deactivate
		'''

rule feature_extraction_DEL:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/filter/{sample}.MAF0.DEL.chr.bed",
		file2="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/{sample}.bam"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.DEL.features"
	params:
		outdir="{outpath}/{sample}/02_map/bqsr/sort.rmdup.bqsr/",
		umap=config['umap'],
		ref=config['reference']
	log:
		"{outpath}/{sample}/03_variants/mutect2/log/{sample}.DEL.features.log"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh && conda activate py3.7.1 
		export PYTHONWARNINGS="ignore" 
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py {input.file1} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1
		conda deactivate
		'''

rule g1000_avail_acess_filter_SNV:
	input:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.SNV.features"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.SNV.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	shell:
		'''
		cat <(head -n 1 {input}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input} | cut -f 2-3) {input}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_INS:
	input:
		feature="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.INS.features"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.INS.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	shell:
		'''
		cat <(head -n 1 {input.feature}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_DEL:
	input:
		feature="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.DEL.features"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.DEL.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	shell:
		'''
		cat <(head -n 1 {input.feature}) <(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' <(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) <(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		'''

rule Prediction_SNV:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.SNV.no1000g.features"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.predictions"
	params:
		refine_beta=config['refine_beta']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Refine {output}
		'''

rule Prediction_INS:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.no1000g.INS.features"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Refine {output}
		'''

rule Prediction_DEL:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.no1000g.DEL.features"
	output:
		"{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Refine {output}
		'''

rule Annotation_SNV_Mosaic:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.predictions"
	output:
		inputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.{ref_version}.input.predictions.txt",
		outputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/{sample}/03_variants/mutect2/feature",
		outputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV"
	shell:
		'''
		grep "mosaic" {input.file1} | cut -f 1,35 | awk 'BEGIN {{ FS = "~" }} ; {{ print $2"\\t"$3"\\t"$3"\\t"$4"\\t"$5"\\t"$6 }}' | tail -n +2 | awk '$6=="mosaic" {{print $0}}' > {output.inputanno}
		cd {params.outdir}
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl {output.inputanno} /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} -buildver {ref_version} -out {params.outputanno} -remove -protocol refGene,dbnsfp42a -operation g,f -nastring . 
		'''

rule Annotation_INS_Mosaic:
	input:
		file1="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.predictions"
	output:
		inputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.{ref_version}.input.predictions.txt",
		mosaic="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.mosaic.predictions.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/{sample}/03_variants/mutect2/feature",
		outputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.ouput.predictions"
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
		file1="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.predictions"
	output:
		inputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.{ref_version}.input.predictions.txt",
		mosaic="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.mosaic.predictions.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/{sample}/03_variants/mutect2/feature",
		outputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.ouput.predictions"
	shell:
		'''
		grep "hap=3" {input.file1} > {output.mosaic}
		cat {output.mosaic} | awk 'BEGIN {{ FS = "~" }} ; {{ print $2"\\t"$3"\\t"$3+length($4)-1"\\t"$4"\\t"$5 }}' | tail -n +2 > {output.inputanno}
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl {output.inputanno} /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} -buildver {ref_version} -out {params.outputanno} -remove -protocol refGene,dbnsfp42a -operation g,f -nastring . 
		'''

rule extrac_SNVsubvcf:
	input:
		vcf="{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.filter.vcf.gz",
		inputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.sub.vcf",
		tmp="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.sub.vcf.tmp",
		bed="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.subvcf.bed",
		dpaf="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.Geno.DP.AF.txt"
	shell:
		'''
		awk '$6=="mosaic"{{print $0}}' {input.inputanno} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2}}' | sort -k1,1 -k2,2n > {output.bed}
		tabix {input.vcf} -R {output.bed} > {output.tmp}
		cat <(zcat {input.vcf} | grep "^#") <(cat {output.tmp}) > {output.vcf}
		cut -f 1,2,4,5,10 {output.tmp} | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		'''

rule extrac_INSsubvcf:
	input:
		vcf="{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.filter.vcf.gz",
		inputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.sub.vcf",
		tmp="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.sub.vcf.tmp",
		bed="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.subvcf.bed",
		dpaf="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.INS.Geno.DP.AF.txt"
	shell:
		'''
		cat {input.inputanno} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2}}' | sort -k1,1 -k2,2n > {output.bed}
		tabix {input.vcf} -R {output.bed} > {output.tmp}
		cat <(zcat {input.vcf} | grep "^#") <(cat {output.tmp}) > {output.vcf}
		cut -f 1,2,4,5,10 {output.tmp} | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		'''
		
rule extrac_DELsubvcf:
	input:
		vcf="{outpath}/{sample}/03_variants/mutect2/result/{sample}.mt2pon.filter.vcf.gz",
		inputanno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.sub.vcf",
		tmp="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.sub.vcf.tmp",
		bed="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.subvcf.bed",
		dpaf="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.DEL.Geno.DP.AF.txt"
	shell:
		'''
		cat {input.inputanno} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2}}' | sort -k1,1 -k2,2n > {output.bed}
		tabix {input.vcf} -R {output.bed} > {output.tmp}
		cat <(zcat {input.vcf} | grep "^#") <(cat {output.tmp}) > {output.vcf}
		cut -f 1,2,4,5,10 {output.tmp} | sed "s/:/\t/g" | cut -f 1-8 | sed 's/,/\t/g' | awk 'BEGIN{{print "Chr\tpos\tREF\tALT\tGenotype\tDPT_REF\tDPT_ALT\tAF\tDPT"}}1' > {output.dpaf}
		'''

rule gnomAD_filter_vcf:
	input:
		vcf="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.sub.vcf",
		gnomad_2_1_1_0_001="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_gnomad211_genome_exome/{ref_version}_gnomad211_exome_genome_afpop_gt_0.001_chr.txt",
		gnomad_2_1_1_0="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_gnomad211_genome_exome/{ref_version}_gnomad211_exome_genome_afpop_gt_0_chr.txt"
	output:
		vcf4="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.sub.MFv4.vcf",
		vcf3="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.sub.MFv3.vcf"
	shell:
		'''
		awk "NR==FNR{{c[$1,$2,$4,$5]++;next}} !(c[$1,$2,$4,$5]) {{print $0}}" {input.gnomad_2_1_1_0_001} {input.vcf} > {output.vcf4}
		awk "NR==FNR{{c[$1,$2,$4,$5]++;next}} !(c[$1,$2,$4,$5]) {{print $0}}" {input.gnomad_2_1_1_0} {input.vcf} > {output.vcf3}
		'''

rule gnomAD_filter_anno:
	input:
		anno="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		vcf4="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.{ref_version}.sub.MFv4.vcf",
		vcf3="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.SNV.{ref_version}.sub.MFv3.vcf"
	output:
		anno4="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV.MFv4.{ref_version}_multianno.txt",
		anno3="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV.MFv3.{ref_version}_multianno.txt"
	shell:
		'''
		cat <(head -n 1 {input.anno}) <(awk 'NR==FNR{{c[$1,$2,$4,$5]++;next}} c[$1,$2,$4,$5]>0{{print $0}}' {input.vcf4} {input.anno}) > {output.anno4}
		cat <(head -n 1 {input.anno}) <(awk 'NR==FNR{{c[$1,$2,$4,$5]++;next}} c[$1,$2,$4,$5]>0{{print $0}}' {input.vcf3} {input.anno}) > {output.anno3}
		'''

rule mark_tier123:
	input:
		anno5="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		anno4="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV.MFv4.{ref_version}_multianno.txt",
		anno3="{outpath}/{sample}/03_variants/mutect2/feature/{sample}.{cov}.annoSNV.MFv3.{ref_version}_multianno.txt"
	output:
		tier="{outpath}/{sample}/03_variants/mutect2/tier/{sample}.{cov}.SNV.tier.{ref_version}_multianno.txt"
	shell:
		'''
		awk 'BEGIN{{OFS="\t"}}NR==FNR{{c[$1,$2,$4,$5]=$0}} NR!=FNR{{if (c[$2,$3,$5,$6]) {{print "tier2\t"$0}} else {{print "tier3\t"$0}}}}' {input.anno4} <(awk 'BEGIN{{OFS="\t"}}NR==FNR{{c[$1,$2,$4,$5]=$0}} NR!=FNR{{if (c[$1,$2,$4,$5]) {{print "tier1\t"$0}} else {{print "tier3\t"$0}}}}' {input.anno3} {input.anno5}) | sed 's/tier2\ttier1/tier1/g' | sed 's/tier2\ttier3/tier2/g' | sed 's/tier3\ttier3/tier3/g' |sed 's/tier3\ttier1/tier1/g' | sed 's/tier1\tChr/Tier\tChr/' > {output.tier}
		'''

rule summary:
	input:
		expand(["{outpath}/{u.sample}/03_variants/mutect2/tier/{u.sample}.{cov}.SNV.tier.{ref_version}_multianno.txt"], outpath=outpath, u=units.itertuples(), cov=cov, ref_version=ref_version)
	output:
		stat="{outpath}/mosaic_summary.txt"
	shell:
		'''
		wc -l {input} > {output.stat}
		'''