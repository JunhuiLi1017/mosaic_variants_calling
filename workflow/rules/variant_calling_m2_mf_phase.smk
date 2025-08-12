rule Prediction_SNV:
	input:
		file1="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/03_feature/{sample}.SNV.no1000g.features"
	output:
		"{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.predictions"
	params:
		phase_model=config['phase_model']
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.phase_model} Phase {output}
		'''

rule Prediction_INS:
	input:
		file1="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/03_feature/{sample}.INS.no1000g.features"
	output:
		"{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Phase {output}
		'''

rule Prediction_DEL:
	input:
		file1="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/03_feature/{sample}.DEL.no1000g.features"
	output:
		"{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.predictions"
	params:
		refine_beta="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/models_trained/insertions_250x.RF.rds"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		'''
		Rscript /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/mosaicforecast/MosaicForecast-master/Prediction.R {input.file1} {params.refine_beta} Phase {output}
		'''

rule extract_bed_snv:
	input:
		prediction_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.predictions"
	output:
		bed_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.{ref_version}.bed.temp_file"
	shell:
		'''
		#awk '$35=="mosaic"{{print $0}}' {input.prediction_snv} > {params.temp_file}
		awk '$35=="hap=3"{{print $0}}' {input.prediction_snv} > {params.temp_file}
		if [ -s {params.temp_file} ]; then
			cat {params.temp_file} | cut -f 1 | cut -d"~" -f2-5 | sed 's/~/\t/g' | grep -v "id" | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2,$3,$4}}' | sort -k1,1 -k2,2n > {output.bed_snv}
		else
			: > {output.bed_snv}
		fi
		rm {params.temp_file}
		'''

rule extract_bed_ins:
	input:
		prediction_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.predictions"
	output:
		bed_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.{ref_version}.bed.temp_file"
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
		prediction_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.predictions"
	output:
		bed_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.{ref_version}.bed"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		temp_file="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.{ref_version}.bed.temp_file"
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

rule extrac_SNVsubvcf:
	input:
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz",
		bed="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.{ref_version}.bed"
	output:
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
		dpaf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.{ref_version}.Geno.DP.AF.txt"
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
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz",
		bed="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.{ref_version}.bed"
	output:
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.{ref_version}.vcf.gz",
		dpaf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.{ref_version}.Geno.DP.AF.txt"
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
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/01_somatic/01_result/{sample}.mt2pon.filter.vcf.gz",
		bed="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.{ref_version}.bed"
	output:
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
		dpaf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.{ref_version}.Geno.DP.AF.txt"
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
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
	output:
		vcf_anno_dbsnp="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.SNV.anno.dbsnp.{ref_version}.vcf.gz",
		vcf_anno_exome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.SNV.anno.exome.{ref_version}.vcf.gz",
		vcf_anno_genome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_anno_genome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz.tbi"
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
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.INS.{ref_version}.vcf.gz",
	output:
		vcf_anno_dbsnp="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.INS.anno.dbsnp.{ref_version}.vcf.gz",
		vcf_anno_exome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.INS.anno.exome.{ref_version}.vcf.gz",
		vcf_anno_genome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_anno_genome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz.tbi"
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
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/04_mosaic_phase/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
	output:
		vcf_anno_dbsnp="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.DEL.anno.dbsnp.{ref_version}.vcf.gz",
		vcf_anno_exome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.DEL.anno.exome.{ref_version}.vcf.gz",
		vcf_anno_genome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_anno_genome="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz.tbi"
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
		vcf_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.SNV.anno.exome.genome.{ref_version}.vcf.gz.tbi",
		vcf_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.INS.anno.exome.genome.{ref_version}.vcf.gz.tbi",
		vcf_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz",
		tbi_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/05_anno_gnomAD_dbsnp_phase/{sample}.{cov}.DEL.anno.exome.genome.{ref_version}.vcf.gz.tbi"
	output:
		vcf_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
		tbi_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.SNV.{ref_version}.vcf.gz.tbi",
		vcf_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.INS.{ref_version}.vcf.gz",
		tbi_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.INS.{ref_version}.vcf.gz.tbi",
		vcf_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
		tbi_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.DEL.{ref_version}.vcf.gz.tbi"
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
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.SNV.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.SNV.{ref_version}.vcf.gz.tbi"
	output:
		outputanno="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt"
	params:
		annovar_dir=config['annovar_dir'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		outputanno="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoSNV"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		if [ $(zgrep -v '^#' {input.vcf} | wc -l) != 0 ]; then
			perl {params.annovar_dir}/table_annovar.pl \
			{input.vcf} \
			{params.annovar_dir}/humandb_{REF_version} \
			-buildver {REF_version} \
			-out {params.outputanno} \
			-remove \
			-protocol refGene,dbnsfp42a,clinvar_20240917 \
			-operation g,f,f \
			-nastring . \
			-vcfinput
		else
		 # VCF is empty, create empty output file
			touch {output.outputanno}
			touch {params.outputanno}.{REF_version}_multianno.vcf
		fi
		'''

rule Annotation_INS:
	input:
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.INS.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.INS.{ref_version}.vcf.gz.tbi"
	output:
		outputanno="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoINS.{ref_version}_multianno.txt"
	params:
		annovar_dir=config['annovar_dir'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		outputanno="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoINS"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		'''
		if [ $(zgrep -v '^#' {input.vcf} | wc -l) != 0 ]; then
			perl {params.annovar_dir}/table_annovar.pl \
			{input.vcf} \
			{params.annovar_dir}/humandb_{REF_version} \
			-buildver {REF_version} \
			-out {params.outputanno} \
			-remove \
			-protocol refGene,dbnsfp42a,clinvar_20240917 \
			-operation g,f,f \
			-nastring . \
			-vcfinput
		else
		 # VCF is empty, create empty output file
			touch {output.outputanno}
			touch {params.outputanno}.{REF_version}_multianno.vcf
		fi
		'''

rule Annotation_DEL:
	input:
		vcf="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.DEL.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/06_mosaic_tier_phase/{sample}.{cov}.DEL.{ref_version}.vcf.gz.tbi"
	output:
		outputanno="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
	params:
		annovar_dir=config['annovar_dir'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000),
		outputanno="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoDEL"
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
			{params.annovar_dir}/humandb_{REF_version} \
			-buildver {REF_version} \
			-out {params.outputanno} \
			-remove \
			-protocol refGene,dbnsfp42a,clinvar_20240917 \
			-operation g,f,f \
			-nastring . \
			-vcfinput
		else
		 # VCF is empty, create empty output file
			touch {output.outputanno}
			touch {params.outputanno}.{REF_version}_multianno.vcf
		fi
		'''

rule reformat_annotation_clinvar:
	input:
		anno_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		anno_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoINS.{ref_version}_multianno.txt",
		anno_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
	output:
		anno_snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		anno_ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.annoINS.{ref_version}_multianno.txt",
		anno_del="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
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
		tier_anno="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt"
	output:
		txt="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annocnv_gnomadlof.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{cov}.{ref_version}.m2.rcnv_lof.log"
	params:
		ref_version=config['ref_version'],
		gnomad_LoF=config['gnomad_LoF'],
		rCNV_gene_score=config['rCNV_gene_score'],
		sub="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.SNV.{ref_version}.exonic_splicing_multianno.txt",
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
		txt="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/07_annovar_phase/{sample}.{cov}.annocnv_gnomadlof.{ref_version}_multianno.txt"
	output:
		txt="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.exonic.splicing.{ref_version}.txt"
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
		snv_exonic="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.exonic.splicing.{ref_version}.txt",
		snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.annoSNV.{ref_version}_multianno.txt",
		ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.annoINS.{ref_version}_multianno.txt",
		deletion="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/08_annovar_reformat_phase/{sample}.{cov}.annoDEL.{ref_version}_multianno.txt"
	output:
		snv="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/{sample}.{cov}.SNV.{ref_version}.mosaic_summary_phase.txt",
		ins="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/{sample}.{cov}.INS.{ref_version}.mosaic_summary_phase.txt",
		deletion="{outpath}/03_variants_case_control/01_mutect2_pon/02_mosaicforecast/{sample}.{cov}.DEL.{ref_version}.mosaic_summary_phase.txt"
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