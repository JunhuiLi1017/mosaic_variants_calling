import os
import argparse

def split_by_chromosome(input_file, output_file, chromosome):
    """
    Filter a tab-delimited file for a specific chromosome, writing to a single output file with header.
    
    Args:
        input_file (str): Path to the input tab-delimited file.
        output_file (str): Path to the output file.
        chromosome (str): Chromosome to filter (e.g., 'chr1').
    """
    # Create output directory if it exists and is non-empty
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Flag to track if any variants are found
    has_variants = False
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read and write the header
        header = infile.readline().strip()
        outfile.write(header + '\n')
        
        # Process each data row
        for line in infile:
            # Skip empty lines
            if not line.strip():
                continue
            
            # Get chromosome from the first column
            chr_name = line.split('\t')[0].strip()
            if not chr_name:
                continue  # Skip rows with empty chromosome
            
            # Write row if it matches the specified chromosome
            if chr_name == chromosome:
                outfile.write(line)
                has_variants = True
    
    # Log result
    if has_variants:
        print(f"Created: {output_file}")
    else:
        print(f"No variants for {chromosome}, created empty file with header: {output_file}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Filter a tab-delimited file for a specific chromosome with header.")
    parser.add_argument("input_file", help="Path to the input tab-delimited file")
    parser.add_argument("-o", "--output-file", required=True, help="Path to the output file")
    parser.add_argument("--chr", required=True, help="Chromosome to filter (e.g., chr1)")
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.isfile(args.input_file):
        parser.error(f"Input file does not exist: {args.input_file}")
    
    # Run the splitting function
    split_by_chromosome(args.input_file, args.output_file, args.chr)

if __name__ == "__main__":
    main()
