rule map_reads:
	input:
		get_fastqc_fastp
	output:
		raw_bam=temp("{outpath}/02_map/01_raw/{sample}.{library}.{flowlane}.raw.bam")
	log:
		"{outpath}/logs/bwa/{sample}.{library}.{flowlane}.bwa.log"
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	params:
		rg=get_reads_group,
		ref=config['reference']
	shell:
		"""
		bwa mem -t {threads} \
		-K 500000000 \
		-M {params.rg} \
		{params.ref} \
		{input[0]} \
		{input[1]} \
		| samtools view -@ {threads} -b -o {output.raw_bam} > {log} 2>&1
		""" 

rule map_reads_sort:
	input:
		"{outpath}/02_map/01_raw/{sample}.{library}.{flowlane}.raw.bam"
	output:
		o2=temp("{outpath}/02_sort/{sample}.{library}.{flowlane}.sort.bam"),
		o3="{outpath}/02_sort/{sample}.{library}.{flowlane}.sort.stat"
	log:
		"{outpath}/logs/bwa/{sample}.{library}.{flowlane}.bwa.sort.log"
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	params:
		tmpdir="{outpath}/02_sort/tmpdir_{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
	container:
		config["sambamba_1.0.1"]
	shell:
		"""
		mkdir -p {params.tmpdir} && \
		sambamba sort \
		-t {threads} \
		-m {params.command_mem}M \
		-o {output.o2} \
		--tmpdir {params.tmpdir} \
		{input} > {log} 2>&1
		samtools stats {output.o2} > {output.o3}
		samtools index {output.o2}
		"""

rule map_reads_sort_index:
	input:
		"{outpath}/02_sort/{sample}.{library}.{flowlane}.sort.bam"
	output:
		"{outpath}/02_sort/{sample}.{library}.{flowlane}.sort.stat"
	log:
		"{outpath}/logs/bwa/{sample}.{library}.{flowlane}.bwa.sort.indelx.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
	container:
		config["samtools_1.20"]
	shell:
		"""
		samtools index {input}
		samtools stats {input} > {output}
		"""

rule map_reads_frag:
	input:
		i1="{outpath}/02_sort/{sample}.{library}.{flowlane}.sort.bam"
	output:
		o1="{outpath}/02_sort/{sample}.{library}.{flowlane}.sort.insert.png",
		o2="{outpath}/02_sort/{sample}.{library}.{flowlane}.sort.fragment.txt"
	log:
		"{outpath}/logs/bwa/{sample}.{library}.{flowlane}.bwa.sort.frag.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["terra_deeptools_0.1"]
	shell:
		"""
		bamPEFragmentSize -b {input.i1} \
		-o {output.o1} \
		--maxFragmentLength 2000 \
		--table {output.o2} > {log} 2>&1
		""" 

rule remove_dup:
	input:
		lambda wildcards: [
			f"{wildcards.outpath}/02_sort/{wildcards.sample}.{u.library}.{u.flowlane}.sort.bam"
			for u in units[units["sample"] == wildcards.sample].itertuples()
		]
	output:
		bam=temp("{outpath}/02_map/03_dup/{sample}.sort.rmdup.bam"),
		bai="{outpath}/02_map/03_dup/{sample}.sort.rmdup.bai",
		mtx="{outpath}/02_map/03_dup/{sample}.sort.rmdup.matrix"
	log:
		"{outpath}/logs/bwa/{sample}.removeduplicate.log"
	params:
		dedup=lambda wildcards: (
			f"--REMOVE_DUPLICATES true" 
			if units[units['sample'] == wildcards.sample]['pcr_based'].any() 
			else "--REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
		),
		tmpdir="{outpath}/02_map/03_dup/tmpdir_{sample}",
		input_args=lambda wildcards, input: " ".join(f"--INPUT {bam}" for bam in input),
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 8000)
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
		#mem_mb=5000
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		mkdir -p {params.tmpdir}
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		MarkDuplicates {params.input_args} \
		--METRICS_FILE {output.mtx} \
		--OUTPUT {output.bam} \
		{params.dedup} \
		--TMP_DIR {params.tmpdir} \
		--CREATE_INDEX true > {log} 2>&1
		"""

rule SetNmMdAndUqTags:
	input:
		bam="{outpath}/02_map/03_dup/{sample}.sort.rmdup.bam"
	output:
		bam=temp("{outpath}/02_map/04_settags/{sample}.sort.rmdup.settags.bam")
	log:
		"{outpath}/02_map/logs/{sample}.settags.log"
	params:
		gatk=config['gatk_current_using'],
		ref=config['reference'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		SetNmMdAndUqTags \
		-I {input.bam} \
		-O {output.bam} \
		-R {params.ref} \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE true  > {log} 2>&1
		"""

rule SetNmMdAndUqTags_sortbam:
	input:
		"{outpath}/02_map/04_settags/{sample}.sort.rmdup.settags.bam"
	output:
		bam="{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam"
	log:
		"{outpath}/02_map/logs/{sample}.sort.log"
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	params:
		tmpdir="{outpath}/02_map/04_settags/tmpdir_{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
	container:
		config["sambamba_1.0.1"]
	shell:
		"""
		mkdir -p {params.tmpdir} && \
		sambamba sort \
		-t {threads} \
		-m {params.command_mem}M \
		-o {output.bam} \
		--tmpdir {params.tmpdir} \
		{input} > {log} 2>&1
		""" 

rule SetNmMdAndUqTags_sortbam_index:
	input:
		"{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam"
	output:
		"{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam.bai"
	log:
		"{outpath}/02_map/logs/{sample}.sort.index.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	params:
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
	container:
		config["samtools_1.20"]
	shell:
		"""
		samtools index {input}
		"""

rule BaseRecalibrator:
	input:
		bam="{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam"
	output:
		o1="{outpath}/02_map/05_BaseRecalibrator_subchr/{sample}.{chr}.recal_data.table"
	log:
		"{outpath}/02_map/logs/{sample}.{chr}.recal.log"
	params:
		gatk=config['gatk_current_using'],
		ref=config['reference'],
		dbsnp138=config['dbsnp138'],
		g1000_known_indels=config['g1000_known_indels'],
		mills_and_1000g=config['mills_and_1000g'],
		interval_list=intervals_dir + "/{chr}.intervals.list",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		BaseRecalibrator \
		-I {input.bam} \
		-O {output.o1} \
		-R {params.ref} \
		-L {params.interval_list} \
		--known-sites {params.dbsnp138} \
		--known-sites {params.g1000_known_indels} \
		--known-sites {params.mills_and_1000g} > {log} 2>&1
		"""

rule GatherBQSRreports:
	input:
		lambda wildcards: [f"{wildcards.outpath}/02_map/05_BaseRecalibrator_subchr/{wildcards.sample}.{chr}.recal_data.table" for chr in CHROMOSOMES]
	output:
		report = "{outpath}/02_map/06_BaseRecalibrator/{sample}.recal_data.table" 
	log:
		"{outpath}/02_map/logs/{sample}.GatherBQSRreports.log"
	params:
		gatk=config['gatk_current_using'],
		input_args=lambda wildcards, input: " ".join(f"-I {report}" for report in input),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		GatherBQSRReports \
		{params.input_args} \
		-O {output.report} > {log} 2>&1
		"""

rule ApplyBQSR:
	input:
		bam="{outpath}/02_map/04_settags/{sample}.rmdup.settags.sort.bam",
		recal_table="{outpath}/02_map/06_BaseRecalibrator/{sample}.recal_data.table"
	output:
		bam="{outpath}/02_map/07_bqsr_subchr/{sample}.{chr}.bam",
		bai="{outpath}/02_map/07_bqsr_subchr/{sample}.{chr}.bai"
	log:
		"{outpath}/02_map/logs/{sample}.{chr}.applyBQSR.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval_list=intervals_dir + "/{chr}.intervals.list",
		prefix="{outpath}/02_map/07_bqsr_subchr/{sample}.{chr}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		ApplyBQSR \
		-I {input.bam} \
		-O {output.bam} \
		-R {params.ref} \
		-L {params.interval_list} \
		--bqsr-recal-file {input.recal_table} > {log} 2>&1
		"""

rule samtools_stats_chr_ApplyBQSR:
	input:
		bam="{outpath}/02_map/07_bqsr_subchr/{sample}.{chr}.bam",
		bai="{outpath}/02_map/07_bqsr_subchr/{sample}.{chr}.bai"
	output:
		stats="{outpath}/02_map/07_bqsr_subchr/{sample}.{chr}.samtools_stats.txt"
	log:
		"{outpath}/02_map/logs/{sample}.{chr}.samtools_stats.log"
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


rule GatherBQSRBam:
	input:
		bam=lambda wildcards: [f"{wildcards.outpath}/02_map/07_bqsr_subchr/{wildcards.sample}.{chr}.bam" for chr in CHROMOSOMES],
		bai=lambda wildcards: [f"{wildcards.outpath}/02_map/07_bqsr_subchr/{wildcards.sample}.{chr}.bai" for chr in CHROMOSOMES]
	output:
		bam="{outpath}/02_map/08_bqsr/{sample}.bam",
		bai="{outpath}/02_map/08_bqsr/{sample}.bai"
	log:
		"{outpath}/02_map/logs/{sample}.gatherbqsrbam.log"
	params:
		gatk=config['gatk_current_using'],
		input_args=lambda wildcards, input: " ".join(f"--INPUT {bam}" for bam in input.bam),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		GatherBamFiles \
		{params.input_args} \
		--OUTPUT {output.bam} \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE true > {log} 2>&1
		"""