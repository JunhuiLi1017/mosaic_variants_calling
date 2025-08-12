import argparse
import re

def parse_genome_regions(region_file):
    """Read the genome region file and return a dictionary of chromosome lengths."""
    chrom_lengths = {}
    with open(region_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Parse format like chr1:1-1002028
            match = re.match(r'(chr[\w\d]+):(\d+)-(\d+)', line)
            if not match:
                raise ValueError(f"Invalid region format in {region_file}: {line}. Expected 'chrX:start-end'.")
            chrom, start, end = match.groups()
            length = int(end) - int(start) + 1
            if length <= 0:
                raise ValueError(f"Invalid region length for {chrom} in {region_file}: start={start}, end={end}")
            chrom_lengths[chrom] = length
    return chrom_lengths

def calculate_mean_depth(input_cov, input_genome_region, output):
    """Calculate mean depth for the whole genome and each chromosome."""
    # Get chromosome lengths
    chrom_lengths = parse_genome_regions(input_genome_region)
    total_genome_length = sum(chrom_lengths.values())
    
    # Initialize per-chromosome coverage and base counts
    chrom_coverage = {chrom: 0 for chrom in chrom_lengths}
    chrom_bases = {chrom: 0 for chrom in chrom_lengths}
    
    # Read COV data
    with open(input_cov, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            chrom, depth, bases = parts
            if chrom not in chrom_lengths:
                print(f"Warning: Chromosome {chrom} in COV data not found in genome regions. Skipping.")
                continue
            try:
                depth = int(depth)
                bases = int(bases)
            except ValueError:
                print(f"Warning: Invalid depth or bases in line: {line}. Skipping.")
                continue
            chrom_coverage[chrom] += depth * bases
            chrom_bases[chrom] += bases
    
    # Calculate per-chromosome mean depths
    chrom_mean_depths = {}
    total_coverage = 0
    total_bases = 0
    for chrom in chrom_lengths:
        zero_coverage_bases = chrom_lengths[chrom] - chrom_bases.get(chrom, 0)
        if zero_coverage_bases < 0:
            raise ValueError(f"Total bases ({chrom_bases.get(chrom, 0)}) exceeds chromosome length ({chrom_lengths[chrom]}) for {chrom}.")
        total_bases_with_zero = chrom_bases.get(chrom, 0) + zero_coverage_bases
        mean_depth = chrom_coverage.get(chrom, 0) / total_bases_with_zero if total_bases_with_zero > 0 else 0
        chrom_mean_depths[chrom] = mean_depth
        total_coverage += chrom_coverage.get(chrom, 0)
        total_bases += chrom_bases.get(chrom, 0)
    
    # Calculate whole-genome mean depth
    zero_coverage_bases = total_genome_length - total_bases
    if zero_coverage_bases < 0:
        raise ValueError(f"Total bases ({total_bases}) exceeds total genome length ({total_genome_length}).")
    total_bases_with_zero = total_bases + zero_coverage_bases
    whole_genome_mean_depth = total_coverage / total_bases_with_zero if total_bases_with_zero > 0 else 0
    
    # Write results to output file
    with open(output, 'w') as f:
        for chrom in sorted(chrom_mean_depths.keys()):
            f.write(f"{chrom}: {chrom_mean_depths[chrom]:.2f}\n")
        f.write(f"all: {whole_genome_mean_depth:.2f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate mean depth from per-chromosome samtools stats COV data and genome regions.")
    parser.add_argument("--input_cov", required=True, help="Path to combined samtools stats COV file (format: chrom depth bases)")
    parser.add_argument("--input_genome_region", required=True, help="Path to file containing chromosome lengths (format: chr length or chr start end)")
    parser.add_argument("--output", required=True, help="Path to output file for mean depth (format: chr1: mean_depth, ..., all: mean_depth)")
    args = parser.parse_args()
    
    calculate_mean_depth(args.input_cov, args.input_genome_region, args.output)
