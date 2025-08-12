rule collect_cov:
	input:
		expand("{{outpath}}/02_map/07_bqsr_subchr/{{sample}}.{chr}.samtools_stats.txt",chr=CHROMOSOMES)
	output:
		"{outpath}/02_map/06_stats/{sample}.cov.txt"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		grep "^COV" {input} | awk -F'/' '{{print $NF}}' | sed "s/\.txt.COV//g" | awk -F'.' '{{print $NF}}' | cut -f 1,3- > {output}
		"""

rule calculate_mean_depth:
	input:
		"{outpath}/02_map/06_stats/{sample}.cov.txt"
	output:
		"{outpath}/02_map/06_stats/{sample}.mean_depth.txt"
	params:
		calculate_mean_depth_script=config["calculate_mean_depth_script"],
		genome_region=intervals_dir + "/all.chr.txt"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["py_3.21"]
	shell:
		"""
		grep "^COV" {input} | awk -F'/' '{{print $NF}}' | sed "s/\.txt.COV//g" | awk -F'.' '{{print $NF}}' | cut -f 1,3- > {output}
		python {params.calculate_mean_depth_script} --input_cov {input} --input_genome_region {params.genome_region} --output {output}
		"""

rule samtools_stats:
	input:
		bam="{outpath}/02_map/08_bqsr/{sample}.bam",
		bai="{outpath}/02_map/08_bqsr/{sample}.bai"
	output:
		stats="{outpath}/02_map/06_stats/{sample}.samtools_stats.txt"
	log:
		"{outpath}/02_map/logs/{sample}.samtools_stats.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["samtools_1.20"]
	shell:
		"""
		samtools stats {input.bam} > {output.stats}
		"""

rule samtools_idxstats:
	input:
		bam="{outpath}/02_map/08_bqsr/{sample}.bam",
		bai="{outpath}/02_map/08_bqsr/{sample}.bai"
	output:
		idxstats="{outpath}/02_map/06_stats/{sample}.samtools_idxstats.txt"
	log:
		"{outpath}/02_map/logs/{sample}.samtools_idxstats.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["samtools_1.20"]
	shell:
		"""
		samtools idxstats {input.bam} > {output.idxstats}
		"""

rule samtools_flagstat:
	input:
		bam="{outpath}/02_map/08_bqsr/{sample}.bam",
		bai="{outpath}/02_map/08_bqsr/{sample}.bai"
	output:
		flagstat="{outpath}/02_map/06_stats/{sample}.samtools_flagstats.txt"
	log:
		"{outpath}/02_map/logs/{sample}.samtools_flagstat.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["samtools_1.20"]
	shell:
		"""
		samtools flagstat {input.bam} > {output.flagstat}
		"""

rule collect_duplicate_metrics:
	input:
		dup_metrics="{outpath}/02_map/03_dup/{sample}.sort.rmdup.matrix"
	output:
		metrics="{outpath}/02_map/06_stats/{sample}.duplicate_metrics.txt"
	params:
		sample="{sample}"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		# Convert GATK MarkDuplicates metrics to MultiQC-compatible format
		echo "Sample	Library	Reads	Duplicates	Duplication_rate" > {output.metrics}
		if [ -f {input.dup_metrics} ]; then
			# Extract metrics from GATK output
			TOTAL_READS=$(grep -A 1 "READ_PAIRS_EXAMINED" {input.dup_metrics} | tail -1 | cut -f 3)
			DUPLICATE_READS=$(grep -A 1 "READ_PAIR_DUPLICATES" {input.dup_metrics} | tail -1 | cut -f 7)
			DUPLICATE_RATE=$(grep -A 1 "PERCENT_DUPLICATION" {input.dup_metrics} | tail -1 | cut -f 9)
			echo "{params.sample}	$TOTAL_READS	$DUPLICATE_READS	$DUPLICATE_RATE" >> {output.metrics}
		else
			echo "{params.sample}	0	0	0.0" >> {output.metrics}
		fi
		"""

rule multiqc_mapping:
	input:
		# Collect all statistics files for MultiQC
		samtools_stats=lambda wildcards: expand(
			"{outpath}/02_map/06_stats/{sample}.samtools_stats.txt",
			outpath=wildcards.outpath,
			sample=Sample
		),
		samtools_idxstats=lambda wildcards: [
			f"{wildcards.outpath}/02_map/06_stats/{sample}.samtools_idxstats.txt"
			for sample in Sample
		],
		samtools_flagstats=lambda wildcards: expand(
			"{outpath}/02_map/06_stats/{sample}.samtools_flagstats.txt",
			outpath=wildcards.outpath,
			sample=Sample
		),
		duplicate_metrics=lambda wildcards: [
			f"{wildcards.outpath}/02_map/06_stats/{sample}.duplicate_metrics.txt"
			for sample in Sample
		],
		coverage_stats=lambda wildcards: expand(
			"{outpath}/02_map/06_stats/{sample}.mean_depth.txt",
			outpath=wildcards.outpath,
			sample=Sample
		)
	output:
		report_map="{outpath}/02_map/multiqc_report_map.html",
		report_raw="{outpath}/02_map/06_stats/multiqc_report.html",
		data_dir="{outpath}/02_map/06_stats/multiqc_data"
	params:
		outdir="{outpath}/02_map/06_stats",
		title="Mapping Quality Report"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config['multiqc_1.22.3']
	shell:
		"""
		multiqc \
			--outdir {params.outdir} \
			--title "{params.title}" \
			--filename multiqc_report.html \
			--data-dir {output.data_dir} \
			--force \
			{params.outdir}
		cp {output.report_raw} {output.report_map}
		"""