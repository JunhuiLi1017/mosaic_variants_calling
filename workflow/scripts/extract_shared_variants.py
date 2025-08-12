import pysam
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Extract shared variants between tumor and normal samples from a VCF file, with './.' matching any genotype and differing from T166.")
parser.add_argument("--input-vcf", required=True, help="Input VCF file (can be compressed)")
parser.add_argument("--output-vcf", required=True, help="Output VCF file with shared variants")
parser.add_argument("--tumor", required=True, help="Tumor sample name")
parser.add_argument("--normal", required=True, help="Normal sample name")
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
required_samples = [args.tumor, args.normal, "T166"]
missing_samples = [s for s in required_samples if s not in vcf_in.header.samples]
if missing_samples:
    vcf_in.close()
    raise ValueError(f"Missing samples in VCF: {', '.join(missing_samples)}")

# Iterate through VCF records
for record in vcf_in:
    # Get genotypes for tumor, normal, and T166
    tumor_gt = record.samples[args.tumor]["GT"]
    normal_gt = record.samples[args.normal]["GT"]
    t166_gt = record.samples["T166"]["GT"]
    
    # Convert genotypes to string format (e.g., "0/1", "./.")
    tumor_gt_str = "/".join("." if x is None else str(x) for x in tumor_gt)
    normal_gt_str = "/".join("." if x is None else str(x) for x in normal_gt)
    t166_gt_str = "/".join("." if x is None else str(x) for x in t166_gt)
    
    # Check if genotypes match for tumor and normal (with "./." as wildcard) and differ from T166
    if (tumor_gt_str == "./." or normal_gt_str == "./." or tumor_gt_str == normal_gt_str) and (tumor_gt_str != t166_gt_str and normal_gt_str != t166_gt_str):
        # Write the entire record to output VCF
        vcf_out.write(record)

# Close VCF files
vcf_in.close()
vcf_out.close()