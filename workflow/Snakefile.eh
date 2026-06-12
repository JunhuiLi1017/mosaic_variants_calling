include: "rules/00_common.smk"
#include: "rules/01_qc.smk"
#include: "rules/02_map_reads.smk"
#include: "rules/02_map_report.smk"
#include: "rules/03_variant_calling_mutect2.smk"
#include: "rules/03_variant_calling_m2_mosaicforecast_anno.smk"
#include: "rules/03_mutect2_mito.smk"
#include: "rules/03_variant_calling_haplotypecaller.smk"
#include: "rules/03_variant_calling_deepvariant.smk"
#include: "rules/03_variant_calling_deepsomatic.smk"
#include: "rules/03_variant_calling_m2_deepmosaic.smk"
#include: "rules/04_filter_anno_clinvar.smk"
include: "rules/06_expansionhunter.smk"

rule all:
    input:
        expand("{outpath}/05_expansionhunter/02_stranger/{sample}.stranger.vcf", outpath=Outpath, sample=Sample)