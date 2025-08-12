rule fastp:
	input:
		unpack(get_fastq)
	output:
		r1=temp("{outpath}/01_multiqc/fastp/{sample}_{library}_{flowlane}.R1.fastq.gz"),
		r2=temp("{outpath}/01_multiqc/fastp/{sample}_{library}_{flowlane}.R2.fastq.gz"),
		j="{outpath}/01_multiqc/fastp/{sample}_{library}_{flowlane}.json",
		h="{outpath}/01_multiqc/fastp/{sample}_{library}_{flowlane}.html"
	log:
		"{outpath}/logs/fastp/{sample}_{library}_{flowlane}.log"
	params:
		trim_expr=lambda wildcards: (
			(
				f"-f {units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_front1']} "
				if (
					(wildcards.sample, wildcards.library, wildcards.flowlane) in units.index
					and pd.notnull(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq1"])
					and units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq1"] != ""
					and int(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_front1']) >= 0
				) else ""
			)
			+
			(
				f"-t {units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_tail1']} "
				if (
					(wildcards.sample, wildcards.library, wildcards.flowlane) in units.index
					and pd.notnull(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq1"])
					and units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq1"] != ""
					and int(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_tail1']) >= 0
				) else ""
			)
			+
			(
				f"-F {units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_front2']} "
				if (
					(wildcards.sample, wildcards.library, wildcards.flowlane) in units.index
					and pd.notnull(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq2"])
					and units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq2"] != ""
					and int(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_front2']) >= 0
				) else ""
			)
			+
			(
				f"-T {units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_tail2']} "
				if (
					(wildcards.sample, wildcards.library, wildcards.flowlane) in units.index
					and pd.notnull(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq2"])
					and units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "fq2"] != ""
					and int(units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), 'trim_tail2']) >= 0
				) else ""
			)
		),
		in_r2=lambda wildcards, input: f"-I {input.r2}" if is_pair_end else "",
		out_r2=lambda wildcards, output: f"-O {output.r2}" if is_pair_end else ""
	message: "fastp quality control"
	threads: resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config['fastp_0.22.0']
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
		get_fastqc_fastp
	output:
		["{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R1_fastqc.html", "{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R2_fastqc.html"],
		["{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R1_fastqc.zip", "{outpath}/01_multiqc/fastqc/{sample}_{library}_{flowlane}.R2_fastqc.zip"]
	log:
		"{outpath}/logs/fastqc/{sample}_{library}_{flowlane}.log"
	params:
		fastqc_out="{outpath}/01_multiqc/fastqc"
	threads: 
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config['fastqc_0.11.9']
	shell:
		"fastqc -o {params.fastqc_out} {input} > {log} 2>&1"

rule multiqc:
	input:
		get_multiqc_fastp
	output:
		"{outpath}/01_multiqc/multiqc_report.html"
	log:
		"{outpath}/logs/multiqc/multiqc.log"
	params:
		out_multiqc="{outpath}/01_multiqc",
		in_fastqc="{outpath}/01_multiqc/fastqc"
	threads: 
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config['multiqc_1.22.3']
	shell:
		"""
		multiqc -o {params.out_multiqc} {params.in_fastqc} --force
		"""