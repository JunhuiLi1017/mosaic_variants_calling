rule gather_vcf:
	input:
		lambda wildcards: [f"{wildcards.outpath}/03_variants/03_deepsomatic/{wildcards.sample}/{wildcards.sample}.{chr}.{wildcards.ref_version}.output.vcf.gz" for chr in chromosomes]
	output:
		vcf="{outpath}/03_variants/03_deepsomatic/01_result/{sample}/{sample}.{ref_version}.output.vcf.gz"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		vcf-concat -f <(ls {input})  |  bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		'''

rule MT2_initial_filter:
	input:
		vcf="{outpath}/03_variants/03_deepsomatic/01_result/{sample}/{sample}.{ref_version}.output.vcf.gz"
	output:
		"{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.deepsomatic.AF0.02.bed"
	params:
		sample="{sample}",
		af_threshold="0.03" if config["pcr_based"] == True else "0.02"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat <(zcat {input.vcf} | grep -v "^#"| grep PASS | cut -f1,2,4,5,10 | sed "s/:/\\t/g" | sed "s/,/\\t/g" | awk '$10>=0.03 && $10<0.4' | grep -v '0|1' | grep -v '1|0') <(zcat {input.vcf} | grep -v "^#" | grep PASS | cut -f1,2,4,5,10 | sed "s/:/\\t/g" | sed "s/,/\\t/g" | awk '$10>={params.af_threshold} && $10<0.4 ' | grep -E "0\|1|1\|0") | cut -f 1-4,8-10 | awk '{{OFS="\\t";print $1,$2-1,$2,$3,$4,\"{params.sample}\",$5,$6,$7}}' > {output}
		'''

rule repeat_filter:
	input:
		dpsomatic="{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.deepsomatic.AF0.02.bed"
	output:
		"{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.dpsomatic.AF0.02.noSegDup.bed"
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
		subtractBed -a {input.dpsomatic} -b {params.segdup} > {output}
		'''

rule annovar_formatter:
	input:
		"{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.dpsomatic.AF0.02.noSegDup.bed"
	output:
		"{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.dpsomatic.AF0.02.ANNOVAR.list"
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
		file1="{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.dpsomatic.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.MAF0.SNV.chr.bed"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"cat {input.file1}|awk '\''{{OFS=\"\\t\";print $1,$2-1,$2,$4,$5,$6}}'\'' | awk '\''length($4)==1 && length($5)==1'\'' |awk 'BEGIN{{OFS=\"\\t\"}} $1=\"chr\"$1' > {output}"


rule MAF0_extraction_INS:
	input:
		file1="{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.dpsomatic.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.MAF0.INS.chr.bed"
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
		file1="{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.dpsomatic.AF0.02.ANNOVAR.list"
	output:
		"{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.MAF0.DEL.chr.bed"
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

rule ln_s_bam:
	input:
		bam="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam",
		bai="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam.bai"
	output:
		bam="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bam",
		bai="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bai"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		ln -s {input.bam} {output.bam}
		ln -s {input.bai} {output.bai}
		'''

rule feature_extraction_SNV:
	input:
		snv_bed="{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.MAF0.SNV.chr.bed",
		bam="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bam",
		bai="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bai"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.SNV.features"
	params:
		outdir="{outpath}/02_map/settags/{sample}/feature_bam/",
		umap=config['umap'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/deepsomatic/{sample}.SNV.{ref_version}.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate py3.7.1
		export PYTHONWARNINGS="ignore" 
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py {input.snv_bed} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1 
		conda deactivate
		'''

rule feature_extraction_INS:
	message:
		"since snakemake cannot create a conda env from ../envs/py3.7.1.yaml, but we would liek to use --use-conda, so we use conda activate"
	input:
		ins_bed="{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.MAF0.INS.chr.bed",
		bam="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bam",
		bai="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bai"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.INS.features"
	params:
		outdir="{outpath}/02_map/settags/{sample}/feature_bam/",
		umap=config['umap'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/deepsomatic/{sample}.INS.{ref_version}.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate py3.7.1
		export PYTHONWARNINGS="ignore" 
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py {input.ins_bed} {output} {params.outdir} {params.ref} {params.umap} 1 bam > {log} 2>&1
		conda deactivate
		'''

rule feature_extraction_DEL:
	input:
		del_bed="{outpath}/03_variants/03_deepsomatic/02_filter/{sample}/{sample}.{ref_version}.MAF0.DEL.chr.bed",
		bam="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bam",
		bai="{outpath}/02_map/settags/{sample}/feature_bam/{sample}.bai"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.DEL.features"
	params:
		outdir="{outpath}/02_map/settags/{sample}/feature_bam/",
		umap=config['umap'],
		ref=config['reference']
	log:
		"{outpath}/03_variants/logs/deepsomatic/{sample}.DEL.{ref_version}.features.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate py3.7.1
		export PYTHONWARNINGS="ignore" 
		python /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/ReadLevel_Features_extraction.py \
		{input.del_bed} \
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
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.SNV.features"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.SNV.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		cat <(head -n 1 {input}) \
		<(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' \
		<(bedtools intersect -a <(sed 's/~/\t/g' {input} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) \
		<(paste <(sed 's/~/\t/g' {input} | cut -f 2-3) {input}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_INS:
	input:
		feature="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.INS.features"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.INS.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		cat <(head -n 1 {input.feature}) \
		<(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' \
		<(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) \
		<(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		'''

rule g1000_avail_acess_filter_DEL:
	input:
		feature="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.DEL.features"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.DEL.no1000g.features"
	params:
		avail_acess_1000g=config['avail_acess_1000g']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		cat <(head -n 1 {input.feature}) \
		<(awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' \
		<(bedtools intersect -a <(sed 's/~/\t/g' {input.feature} | awk '{{print $2"\t"$3"\t"$3}}' | grep -v "conflict_num") -b {params.avail_acess_1000g} -wa) \
		<(paste <(sed 's/~/\t/g' {input.feature} | cut -f 2-3) {input.feature}) | cut -f 3-) > {output}
		'''

rule Prediction_SNV:
	input:
		file1="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.SNV.no1000g.features"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.{ref_version}.SNV.predictions"
	params:
		refine_beta=config['refine_beta']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R \
		{input.file1} \
		{params.refine_beta} \
		Refine \
		{output}
		'''

rule Prediction_INS:
	input:
		file1="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.INS.no1000g.features"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.{ref_version}.INS.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R \
		{input.file1} \
		{params.refine_beta} \
		Refine \
		{output}
		'''

rule Prediction_DEL:
	input:
		file1="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{ref_version}.DEL.no1000g.features"
	output:
		"{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.{ref_version}.DEL.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R \
		{input.file1} \
		{params.refine_beta} \
		Refine \
		{output}
		'''

rule Annotation_SNV_Mosaic:
	input:
		file1="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.{ref_version}.SNV.predictions"
	output:
		inputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.input.predictions.txt",
		outputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}",
		outputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.annoSNV"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		grep "mosaic" {input.file1} | cut -f 1,35 | awk 'BEGIN {{ FS = "~" }} ; {{ print $2"\\t"$3"\\t"$3"\\t"$4"\\t"$5"\\t"$6 }}' | tail -n +2 | awk '$6=="mosaic" {{print $0}}' > {output.inputanno}
		cd {params.outdir}
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl \
		{output.inputanno} \
		/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} \
		-buildver {ref_version} \
		-out {params.outputanno} \
		-remove -protocol refGene,dbnsfp42a \
		-operation g,f \
		-nastring . 
		'''

rule Annotation_INS_Mosaic:
	input:
		file1="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.{ref_version}.INS.predictions"
	output:
		inputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.{ref_version}.input.predictions.txt",
		mosaic="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.mosaic.predictions.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}",
		outputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.ouput.predictions"
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
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl \
		{output.inputanno} \
		/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} \
		-buildver {ref_version} \
		-out {params.outputanno} \
		-remove -protocol refGene,dbnsfp42a \
		-operation g,f \
		-nastring . 
		'''

rule Annotation_DEL_Mosaic:
	input:
		file1="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.{ref_version}.DEL.predictions"
	output:
		inputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.{ref_version}.input.predictions.txt",
		mosaic="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.mosaic.predictions.{ref_version}_multianno.txt"
	params:
		outdir="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}",
		outputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.ouput.predictions"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		grep "hap=3" {input.file1} > {output.mosaic}
		cat {output.mosaic} | awk 'BEGIN {{ FS = "~" }} ; {{ print $2"\\t"$3"\\t"$3+length($4)-1"\\t"$4"\\t"$5 }}' | tail -n +2 > {output.inputanno}
		perl /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl \
		{output.inputanno} \
		/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_{ref_version} 
		-buildver {ref_version} \
		-out {params.outputanno} \
		-remove -protocol refGene,dbnsfp42a \
		-operation g,f \
		-nastring . 
		'''

rule extrac_SNVsubvcf:
	input:
		vcf="{outpath}/03_variants/03_deepsomatic/01_result/{sample}/{sample}.{ref_version}.output.vcf.gz",
		inputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf",
		tmp="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf.tmp",
		bed="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.subvcf.bed",
		dpaf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.Geno.DP.AF.txt"
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
		vcf="{outpath}/03_variants/03_deepsomatic/01_result/{sample}/{sample}.{ref_version}.output.vcf.gz",
		inputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.{ref_version}.sub.vcf",
		tmp="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.{ref_version}.sub.vcf.tmp",
		bed="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.{ref_version}.subvcf.bed",
		dpaf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.{ref_version}.Geno.DP.AF.txt"
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
		vcf="{outpath}/03_variants/03_deepsomatic/01_result/{sample}/{sample}.{ref_version}.output.vcf.gz",
		inputanno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.{ref_version}.input.predictions.txt"
	output:
		vcf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.{ref_version}.sub.vcf",
		tmp="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.{ref_version}.sub.vcf.tmp",
		bed="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.{ref_version}.subvcf.bed",
		dpaf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.{ref_version}.Geno.DP.AF.txt"
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
		vcf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf",
	output:
		vcf_tier3="{outpath}/03_variants/03_deepsomatic/tier/{sample}/overlap_gnomAD/{sample}.{cov}.SNV.{ref_version}.sub.tier3.{chr}.vcf",
		vcf_tier2="{outpath}/03_variants/03_deepsomatic/tier/{sample}/overlap_gnomAD/{sample}.{cov}.SNV.{ref_version}.sub.tier2.{chr}.vcf"
	params:
		gnomad_2_1_1_0_001="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_gnomad211_genome_exome/sub_chr/{ref_version}_gnomad211_exome_genome_afpop_gt_0.001_chr_snponly_{chr}.txt",
		gnomad_2_1_1_0="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar/humandb_gnomad211_genome_exome/sub_chr/{ref_version}_gnomad211_exome_genome_afpop_gt_0_le_0.001_chr_snponly_{chr}.txt"
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
		vcf_tier3=lambda wildcards: [f"{wildcards.outpath}/03_variants/03_deepsomatic/tier/{wildcards.sample}/overlap_gnomAD/{wildcards.sample}.{wildcards.cov}.SNV.{wildcards.ref_version}.sub.tier3.{chr}.vcf" for chr in chromosomes],
		vcf_tier2=lambda wildcards: [f"{wildcards.outpath}/03_variants/03_deepsomatic/tier/{wildcards.sample}/overlap_gnomAD/{wildcards.sample}.{wildcards.cov}.SNV.{wildcards.ref_version}.sub.tier2.{chr}.vcf" for chr in chromosomes]		
	output:
		vcf_tier3="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier3.vcf",
		vcf_tier2="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier2.vcf"
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
		vcf="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.vcf",
		vcf_tier3="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier3.vcf",
		vcf_tier2="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier2.vcf"
	output:
		vcf_tier1="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier1.vcf"
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
		anno="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		vcf_tier3="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier3.vcf",
		vcf_tier2="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier2.vcf",
		vcf_tier1="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.{ref_version}.sub.tier1.vcf"
	output:
		tier_anno="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.tier.{ref_version}_multianno.txt"
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
		#expand(["{outpath}/03_variants/03_deepsomatic/tier/{u.sample}/{u.sample}.{cov}.SNV.tier.{ref_version}_multianno.txt"], outpath=outpath, u=units.itertuples(), cov=cov, ref_version=ref_version)
		snv="{outpath}/03_variants/03_deepsomatic/tier/{sample}/{sample}.{cov}.SNV.tier.{ref_version}_multianno.txt",
		ins="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.INS.{ref_version}.Geno.DP.AF.txt",
		del1="{outpath}/03_variants/03_deepsomatic/03_feature/{sample}/{sample}.{cov}.DEL.{ref_version}.Geno.DP.AF.txt"
	output:
		snv="{outpath}/03_variants/03_deepsomatic/{sample}.{cov}.SNV.tier.{ref_version}.mosaic_summary.txt",
		ins="{outpath}/03_variants/03_deepsomatic/{sample}.{cov}.INS.{ref_version}.mosaic_summary.txt",
		del1="{outpath}/03_variants/03_deepsomatic/{sample}.{cov}.DEL.{ref_version}.mosaic_summary.txt"
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
