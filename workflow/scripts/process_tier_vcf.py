#!/usr/bin/env python3
import os
import pysam
import sys

def process_tier_vcf(input_vcf, tier_vcf):
    # Open input VCF
    vcf_in = pysam.VariantFile(input_vcf)
    
    # Add TIER header
    vcf_in.header.info.add('TIER', number=1, type='Integer', description='Tier based on AF_genome or AF_exome: 3 (>0.001), 2 (0-0.001], 1 (0 or not found)')
    
    # Determine output mode (compressed or uncompressed)
    output_mode = 'wz' if tier_vcf.endswith('.gz') else 'w'
    
    # Open output VCF
    vcf_out = pysam.VariantFile(tier_vcf, output_mode, header=vcf_in.header)
    
    # Process each variant
    for record in vcf_in:
        # Extract AF_genome and AF_exome
        af_genome = record.info.get('AF_genome', 0.0)
        af_exome = record.info.get('AF_exome', 0.0)
        
        # Convert to float, handle invalid values
        try:
            af_genome = float(af_genome)
        except (ValueError, TypeError):
            af_genome = 0.0
        try:
            af_exome = float(af_exome)
        except (ValueError, TypeError):
            af_exome = 0.0
        
        # Assign tier based on conditions
        if (af_genome > 0.001) or (af_exome > 0.001):
            tier = 3
        elif (0 < af_genome <= 0.001) or (0 < af_exome <= 0.001):
            tier = 2
        else:
            tier = 1
        
        # Add TIER to INFO
        record.info['TIER'] = tier
        
        # Write to output
        vcf_out.write(record)
    
    vcf_in.close()
    vcf_out.close()
    
    # Index compressed output
    if tier_vcf.endswith('.gz'):
        # Use CSI to support large references
        pysam.tabix_index(tier_vcf, preset='vcf', force=True, csi=True)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: annotate_tier.py input.vcf tier.vcf")
        sys.exit(1)
    process_tier_vcf(sys.argv[1], sys.argv[2])
