rule map_reads:
	input:
		get_fastp_fastq
	output:
		raw_bam=temp("{outpath}/02_map/bwa/raw/{sample}/{sample}.{library}.{flowlane}.raw.bam")
	log:
		"{outpath}/logs/bwa/{sample}/{sample}.{library}.{flowlane}.bwa.log"
	threads:
		resource['resource']['huge']['threads']
	resources:
		mem_mb=resource['resource']['huge']['mem_mb']
	params:
		rg=get_reads_group,
		ref=config['reference']
	conda:
		"../envs/bwa.yaml"
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
		"{outpath}/02_map/bwa/raw/{sample}/{sample}.{library}.{flowlane}.raw.bam"
	output:
		o2="{outpath}/02_map/bwa/sort/{sample}/{sample}.{library}.{flowlane}.sort.bam",
		o3="{outpath}/02_map/bwa/sort/{sample}/{sample}.{library}.{flowlane}.sort.stat"
	log:
		"{outpath}/logs/bwa/{sample}/{sample}.{library}.{flowlane}.bwa.sort.log"
	threads:
		resource['resource']['huge']['threads']
	resources:
		mem_mb=resource['resource']['huge']['mem_mb']
	params:
		tmpdir="{outpath}/02_map/bwa/sort/{sample}/tmpdir_{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
	conda:
		"../envs/sambamba1.0.1.yaml"
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

rule map_reads_frag:
	input:
		i1="{outpath}/02_map/bwa/sort/{sample}/{sample}.{library}.{flowlane}.sort.bam"
	output:
		o1="{outpath}/02_map/bwa/sort/{sample}/{sample}.{library}.{flowlane}.sort.insert.png",
		o2="{outpath}/02_map/bwa/sort/{sample}/{sample}.{library}.{flowlane}.sort.fragment.txt"
	log:
		"{outpath}/logs/bwa/{sample}/{sample}.{library}.{flowlane}.bwa.sort.frag.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/deeptools.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate deeptools
		bamPEFragmentSize -b {input.i1} \
		-o {output.o1} \
		--maxFragmentLength 2000 \
		--table {output.o2} > {log} 2>&1
		conda deactivate
		""" 

rule remove_dup:
	input:
		lambda wildcards: [
			f"{outpath}/02_map/bwa/sort/{wildcards.sample}/{wildcards.sample}.{u.library}.{u.flowlane}.sort.bam"
			for u in units[units["sample"] == wildcards.sample].itertuples()
		]
	output:
		bam=temp("{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.bam"),
		bai="{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.bai",
		mtx="{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.matrix"
	log:
		"{outpath}/logs/bwa/{sample}/{sample}.removeduplicate.log"
	params:
		dedup="--REMOVE_DUPLICATES true" if config["pcr_based"] == True else "--REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500",
		tmpdir="{outpath}/02_map/dup/{sample}/tmpdir_{sample}",
		input_args=lambda wildcards, input: " ".join(f"--INPUT {bam}" for bam in input),
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 8000)
	threads:
		resource['resource']['very_high']['threads']
	resources:
		#mem_mb=resource['resource']['very_high']['mem_mb']
		mem_mb=5000
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		mkdir -p {params.tmpdir}
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
		MarkDuplicates {params.input_args} \
		--METRICS_FILE {output.mtx} \
		--OUTPUT {output.bam} \
		{params.dedup} \
		--TMP_DIR {params.tmpdir} \
		--CREATE_INDEX true > {log} 2>&1
		conda deactivate
		"""

rule SetNmMdAndUqTags:
	input:
		bam="{outpath}/02_map/dup/{sample}/{sample}.sort.rmdup.bam"
	output:
		bam="{outpath}/02_map/settags/{sample}/{sample}.sort.rmdup.settags.bam"
	log:
		"{outpath}/02_map/logs/{sample}/{sample}.settags.log"
	params:
		gatk=config['gatk_current_using'],
		ref=config['reference'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} SetNmMdAndUqTags \
		-I {input} \
		-O {output.bam} \
		-R {params.ref} \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE true  > {log} 2>&1
		conda deactivate
		"""

rule sortbam:
	input:
		"{outpath}/02_map/settags/{sample}/{sample}.sort.rmdup.settags.bam"
	output:
		bam="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam",
		bai="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam.bai"
	log:
		"{outpath}/02_map/logs/{sample}/{sample}.sort.log"
	threads:
		resource['resource']['huge']['threads']
	resources:
		mem_mb=resource['resource']['huge']['mem_mb']
	params:
		tmpdir="{outpath}/02_map/settags/{sample}/tmpdir_{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000) // threads
	conda:
		"../envs/sambamba1.0.1.yaml"
	shell:
		"""
		mkdir -p {params.tmpdir} && \
		sambamba sort \
		-t {threads} \
		-m {params.command_mem}M \
		-o {output.bam} \
		--tmpdir {params.tmpdir} \
		{input} > {log} 2>&1
		samtools index {output.bam}
		""" 

rule BaseRecalibrator:
	input:
		bam="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam"
	output:
		o1="{outpath}/02_map/bqsr/{sample}/sub_recal/{sample}.{individual_chr}.recal_data.table"
	log:
		"{outpath}/02_map/logs/{sample}/{sample}.{individual_chr}.recal.log"
	params:
		gatk=config['gatk_current_using'],
		ref=config['reference'],
		dpsnp138=config['dpsnp138'],
		known_indels=config['known_indels'],
		Mills_and_1000G=config['Mills_and_1000G'],
		interval_list=intervals_dir + "/{individual_chr}.intervals.list",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} BaseRecalibrator \
		-I {input} \
		-O {output.o1} \
		-R {params.ref} \
		-L {params.interval_list} \
		--known-sites {params.dpsnp138} \
		--known-sites {params.known_indels} \
		--known-sites {params.Mills_and_1000G} > {log} 2>&1
		conda deactivate
		"""

rule GatherBQSRreports:
	input:
		lambda wildcards: [f"{wildcards.outpath}/02_map/bqsr/{wildcards.sample}/sub_recal/{wildcards.sample}.{chr}.recal_data.table" for chr in chromosomes]
	output:
		report = "{outpath}/02_map/bqsr/{sample}/{sample}.recal_data.table" 
	log:
		"{outpath}/02_map/logs/{sample}/{sample}.GatherBQSRreports.log"
	params:
		gatk=config['gatk_current_using'],
		input_args=lambda wildcards, input: " ".join(f"-I {report}" for report in input),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} GatherBQSRReports \
		{params.input_args} \
		-O {output.report} > {log} 2>&1
		conda deactivate
		"""

rule ApplyBQSR:
	input:
		bam="{outpath}/02_map/settags/{sample}/{sample}.rmdup.settags.sort.bam",
		recal_table="{outpath}/02_map/bqsr/{sample}/{sample}.recal_data.table"
	output:
		bam="{outpath}/02_map/bqsr/{sample}/sub_bqsr/{sample}.{individual_chr}.bam",
		bai="{outpath}/02_map/bqsr/{sample}/sub_bqsr/{sample}.{individual_chr}.bai"
	log:
		"{outpath}/02_map/logs/{sample}.{individual_chr}.applyBQSR.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval_list=intervals_dir + "/{individual_chr}.intervals.list",
		prefix="{outpath}/02_map/bqsr/{sample}/sub_bqsr/{sample}.{individual_chr}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} ApplyBQSR \
		-I {input.bam} \
		-O {output.bam} \
		-R {params.ref} \
		-L {params.interval_list} \
		--bqsr-recal-file {input.recal_table} > {log} 2>&1
		samtools stats {output.bam} > {params.prefix}.txt
		conda deactivate
		"""

rule GatherBQSRBam:
	input:
		bam=lambda wildcards: [f"{wildcards.outpath}/02_map/bqsr/{wildcards.sample}/sub_bqsr/{wildcards.sample}.{chr}.bam" for chr in chromosomes],
		bai=lambda wildcards: [f"{wildcards.outpath}/02_map/bqsr/{wildcards.sample}/sub_bqsr/{wildcards.sample}.{chr}.bai" for chr in chromosomes]
	output:
		bam="{outpath}/02_map/bqsr/{sample}/{sample}.bam",
		bai="{outpath}/02_map/bqsr/{sample}/{sample}.bai"
	log:
		"{outpath}/02_map/logs/{sample}.gatherbqsrbam.log"
	params:
		gatk=config['gatk_current_using'],
		input_args=lambda wildcards, input: " ".join(f"--INPUT {bam}" for bam in input.bam),
		prefix="{outpath}/02_map/bqsr/{sample}/{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} GatherBamFiles \
		{params.input_args} \
		--OUTPUT {output.bam} \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE true > {log} 2>&1
		samtools stats {output.bam} > {params.prefix}.txt
		conda deactivate
		"""