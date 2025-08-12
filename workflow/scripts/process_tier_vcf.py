#!/usr/bin/env python3
import pysam
import sys

def process_tier_vcf(input_vcf, tier_vcf):
    # Open input VCF
    vcf_in = pysam.VariantFile(input_vcf)
    
    # Add TIER header
    vcf_in.header.info.add('TIER', number=1, type='Integer', description='Tier based on AF_genome or AF_exome:4 (>0.05) 3 (0.001-0.05), 2 (0-0.001], 1 (0 or not found)')
    
    # Open output VCF
    vcf_out = pysam.VariantFile(tier_vcf, 'w', header=vcf_in.header)
    
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
        if af_genome > 0.05:
            tier = 4
        elif (0.05 >= af_genome > 0.001) or (0.05 >= af_exome > 0.001):
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

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: annotate_tier.py input.vcf tier.vcf")
        sys.exit(1)
    process_tier_vcf(sys.argv[1], sys.argv[2])
