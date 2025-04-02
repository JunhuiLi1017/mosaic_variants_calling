rule fastp:
	input:
		unpack(get_fastq)
	output:
		r1=temp(f"{outpath}/01_multiqc/fastp/{{sample}}_{{library}}_{{flowlane}}.R1.fastq.gz"),
		r2=temp(f"{outpath}/01_multiqc/fastp/{{sample}}_{{library}}_{{flowlane}}.R2.fastq.gz") if paired_end else [],
		j=f"{outpath}/01_multiqc/fastp/{{sample}}_{{library}}_{{flowlane}}.json",
		h=f"{outpath}/01_multiqc/fastp/{{sample}}_{{library}}_{{flowlane}}.html"
	log:
		f"{outpath}/logs/fastp/{{sample}}_{{library}}_{{flowlane}}.log"
	params:
		trim_expr=f"-f {config['trim_value_f']} -F {config['trim_value_F']}" if config.get('trim_f') and config.get('trim_F') else f"-f {config['trim_value_f']}" if config.get('trim_f') else f"-F {config['trim_value_F']}" if config.get('trim_F') else "",
		in_r2=lambda wildcards, input: f"-I {input.r2}" if paired_end else "",
		out_r2=lambda wildcards, output: f"-O {output.r2}" if paired_end else ""
	message: "fastp quality control"
	threads: resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/fastp.yaml" 
	shell:
		"""
		fastp {params.trim_expr} \
			-i {input.r1} \
			{params.in_r2} \
			-o {output.r1} \
			{params.out_r2} \
			-j {output.j} -h {output.h} \
			--thread {threads} \
			> {log} 2>&1
		"""


rule fastqc:
	input:
		get_fastp_fastq
	output:
		["{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R1_fastqc.html", "{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R2_fastqc.html"] if paired_end else "{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R1_fastqc.html",
		["{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R1_fastqc.zip", "{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R2_fastqc.zip"] if paired_end else "{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R1_fastqc.zip"
	log:
		"{outpath}/logs/fastqc/{sample}_{library}_{flowlane}.log"
	threads: resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/fastqc.yaml" 
	shell:
		"fastqc -o {outpath}/01_multiqc/fastqc {input} > {log} 2>&1"

rule samtools_stats:
	input:
		"{outpath}/02_map/bqsr/{sample}.sort.rmdup.bqsr.bam"
	output:
		"{outpath}/02_map/bqsr_stat/{sample}.sort.rmdup.bqsr.stat"
	log:
		 "{outpath}/logs/bwa/{sample}.bqsr.stat.log"
	shell:
		"samtools stats {input} > {output}"

rule multiqc:
	input:
		get_multiqc_input
	output:
		"{outpath}/01_multiqc/multiqc_report.html"
	conda:
		"../envs/multiqc.yaml"
	log:
		"{outpath}/logs/multiqc/multiqc.log"
	threads: resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		multiqc -o {outpath}/01_multiqc {outpath}/01_multiqc/fastqc --force
		"""