import pysam
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Extract shared or unique variants between tumor and normal samples from a VCF file.")
parser.add_argument("--input-vcf", required=True, help="Input VCF file (can be compressed)")
parser.add_argument("--output-vcf", required=True, help="Output VCF file with variants")
parser.add_argument("--tumor", required=True, help="Tumor sample name")
parser.add_argument("--normal", required=True, help="Normal sample name")
parser.add_argument("--type", choices=["shared", "uniq"], default="shared", help="Type of variants to extract: 'shared' (tumor and normal match, differ from T166) or 'uniq' (tumor and normal differ)")
args = parser.parse_args()

# Open input VCF
try:
    vcf_in = pysam.VariantFile(args.input_vcf, "r")
except Exception as e:
    raise ValueError(f"Failed to open input VCF {args.input_vcf}: {e}")

# Create output VCF with the same header as input
try:
    vcf_out = pysam.VariantFile(args.output_vcf, "w", header=vcf_in.header)
except Exception as e:
    vcf_in.close()
    raise ValueError(f"Failed to create output VCF {args.output_vcf}: {e}")

# Verify that required samples exist in the VCF
required_samples = [args.tumor, args.normal]
if args.type == "shared":
    required_samples.append("T166")
missing_samples = [s for s in required_samples if s not in vcf_in.header.samples]
if missing_samples:
    vcf_in.close()
    raise ValueError(f"Missing samples in VCF: {', '.join(missing_samples)}")

# Iterate through VCF records
for record in vcf_in:
    # Get genotypes for tumor and normal
    tumor_gt = record.samples[args.tumor]["GT"]
    normal_gt = record.samples[args.normal]["GT"]
    
    # Convert genotypes to string format (e.g., "0/1", "./.")
    tumor_gt_str = "/".join("." if x is None else str(x) for x in tumor_gt)
    normal_gt_str = "/".join("." if x is None else str(x) for x in normal_gt)
    
    if args.type == "shared":
        # Get T166 genotype for shared mode
        t166_gt = record.samples["T166"]["GT"]
        t166_gt_str = "/".join("." if x is None else str(x) for x in t166_gt)
        
        # Check if genotypes match for tumor and normal and differ from T166
        if tumor_gt_str == normal_gt_str and tumor_gt_str != t166_gt_str:
            vcf_out.write(record)
    else:  # args.type == "uniq"
        # Check if tumor and normal genotypes differ
        if tumor_gt_str != normal_gt_str:
            vcf_out.write(record)

# Close VCF files
vcf_in.close()
vcf_out.close()
