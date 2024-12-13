#bsub -W 96:00 -q long 'source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake; bash submit.sh &> submit.log'
source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake
snakemake -s /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/mosaic_variants_calling/workflow/Snakefile --configfile /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/mosaic_variants_calling/config/config.yaml -p -j 99 --latency-wait 500 --cluster 'bsub -q long -o /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/mosaic_variants_calling/workflow/main.log -R "rusage[mem=2000]" -n 8 -R span[hosts=1] -W 96:00'
conda deactivate
