#bsub -W 72:00 -q long 'source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake; bash submit.sh &> submit.log'
source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake
snakemake -s /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/mosaic_variants_calling/workflow/Snakefile --configfile config/config.yaml -p -j 99 --latency-wait 500 --default-resources mem_mb=10000 disk_mb=10000 --cluster 'bsub -q long -o main.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 72:00'
conda deactivate
