#bsub -W 72:00 -q long 'source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake; bash submit.sh &> submit.log'
source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake
snakemake -s /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/mosaic_variants_calling/workflow/Snakefile --configfile ../config/config.hg38.yaml -p -j 99 --latency-wait 500 --cluster 'bsub -q long -o main_hg38_NC19-23.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 72:00' --use-singularity
conda deactivate