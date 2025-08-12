import pandas as pd
import re
import argparse

def reformat_vcf_dataset(input_file, output_file, info_column='Otherinfo12', gt_column='Otherinfo13', gt2_column=None, gnomad_column='Otherinfo11'):
    """
    Reformat a VCF-like dataset into a tab-delimited file with specified columns,
    appending all remaining input columns except geneName and Gene. Includes Tier if 
    present in input, otherwise omits it. Extracts GT, AD, DP from gt_column and 
    optionally gt2_column based on field positions in info_column, calculates AF as 
    depth_alt/(depth_alt+depth_ref) when valid. Includes the specified gnomad_column 
    and extracts AF_exome, AF_popmax_exome, AF_genome, AF_popmax_genome if present, 
    otherwise outputs "not_found_in_gnomad". Includes ClinVar columns CLNALLELEID, 
    CLNDN, CLNDISDB, CLNREVSTAT, CLNSIG, ONCDN, ONCDISDB, ONCREVSTAT, ONC, SCIDN, 
    SCIDISDB, SCIREVSTAT, SCI. If pHaplo, pTriplo, _loeuf, _pli columns are present, 
    they are appended after TIER in the output. If gt2_column is specified, includes 
    additional columns AF_2, DP_2, Depth_Ref_2, Depth_Alt_2, Genotype_2 for the second 
    genotype field. Handles empty input files (only headers) by producing an output 
    file with the expected columns but no data rows.
    
    Args:
        input_file (str): Path to the input tab-delimited file.
        output_file (str): Path to the output tab-delimited file.
        info_column (str): Name of the column containing VCF field order (default: 'Otherinfo12').
        gt_column (str): Name of the column containing first VCF genotype field (default: 'Otherinfo13').
        gt2_column (str, optional): Name of the column containing second VCF genotype field (default: None).
        gnomad_column (str): Name of the column containing VCF INFO fields (default: 'Otherinfo11').
    
    Returns:
        None: Writes the reformatted data to output_file.
    """
    # Read the tab-delimited input file
    try:
        df = pd.read_csv(input_file, sep="\t", low_memory=False)
        if len(df) == 0:
            print(f"Warning: Input file {input_file} only contains headers. Creating output with headers only.")
            df = pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"Warning: Input file {input_file} is completely empty. Creating output with headers only.")
        df = pd.DataFrame()

    # Define required columns
    required_columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', gnomad_column, info_column, gt_column]
    if gt2_column is not None:
        required_columns.append(gt2_column)

    # If DataFrame is empty, create output with expected columns
    if df.empty:
        output_columns = [
            'Chr', 'Start', 'End', 'Ref', 'Alt', 'AF', 'DP', 'Depth_Ref', 'Depth_Alt', 'Genotype',
            'rsID', 'AF_exome', 'AF_popmax_exome', 'AF_genome', 'AF_popmax_genome', 'TIER',
            'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'ONCDN', 'ONCDISDB',
            'ONCREVSTAT', 'ONC', 'SCIDN', 'SCIDISDB', 'SCIREVSTAT', 'SCI'
        ]
        if gt2_column:
            output_columns.extend(['AF_2', 'DP_2', 'Depth_Ref_2', 'Depth_Alt_2', 'Genotype_2'])
        # Add optional columns
        optional_output_names = ['pHap', 'pTrio', 'leouf', 'pLI']
        output_columns.extend([col for col in optional_output_names])
        # Create empty DataFrame with columns
        output_df = pd.DataFrame(columns=output_columns)
        output_df.to_csv(output_file, sep="\t", index=False)
        print(f"Reformatted empty data written to {output_file}")
        return

    # Check for missing required columns
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    # Function to extract AF, DP, Depth_Ref, Depth_Alt, and Genotype from VCF fields
    def extract_vcf_info(row, genotype_col):
        depth_ref = 0
        depth_alt = 0
        dp = "."
        af = "."
        gt = "."
        
        info_field = row.get(info_column, '')
        field_indices = {}
        if pd.notna(info_field):
            fields = str(info_field).split(':')
            field_indices = {field: idx for idx, field in enumerate(fields)}

        gt_field = row.get(genotype_col, '')
        if pd.notna(gt_field):
            values = str(gt_field).split(':')
            
            if 'AD' in field_indices:
                ad_idx = field_indices['AD']
                if ad_idx < len(values):
                    ad_match = re.match(r'(\d+),(\d+)', values[ad_idx])
                    if ad_match:
                        depth_ref = int(ad_match.group(1))
                        depth_alt = int(ad_match.group(2))
            
            if 'DP' in field_indices:
                dp_idx = field_indices['DP']
                if dp_idx < len(values):
                    dp = values[dp_idx]
            
            if 'GT' in field_indices:
                gt_idx = field_indices['GT']
                if gt_idx < len(values):
                    gt = values[gt_idx]

            if 'VAF' in field_indices or 'AF' in field_indices:
                af_idx = field_indices.get('VAF', field_indices.get('AF'))
                if af_idx < len(values):
                    af = float(values[af_idx])

        if dp == "." and depth_ref + depth_alt > 0:
            dp = depth_ref + depth_alt
        
        if af == "." and depth_ref + depth_alt > 0:
            af = float(depth_alt) / (depth_ref + depth_alt)
        
        return pd.Series({
            'AF': af,
            'DP': dp,
            'Depth_Ref': depth_ref if depth_ref > 0 else ".",
            'Depth_Alt': depth_alt if depth_alt > 0 else ".",
            'Genotype': gt
        })

    # Function to extract gnomAD fields from gnomad_column
    def extract_gnomad_fields(row):
        gnomad_fields = {
            'AF_exome': '.',
            'AF_popmax_exome': '.',
            'AF_genome': '.',
            'AF_popmax_genome': '.',
            'TIER': '.'
        }
        
        info11 = row.get(gnomad_column, '')
        if pd.notna(info11):
            fields = str(info11).split(';')
            for field in fields:
                if '=' in field:
                    key, value = field.split('=', 1)
                    if key in gnomad_fields:
                        gnomad_fields[key] = value
        
        return pd.Series(gnomad_fields)

    # Initialize output DataFrame with mandatory columns
    output_columns2 = {
        'CLNALLELEID': df.get('CLNALLELEID', pd.Series(['.'] * len(df))),
        'CLNDN': df.get('CLNDN', pd.Series(['.'] * len(df))),
        'CLNDISDB': df.get('CLNDISDB', pd.Series(['.'] * len(df))),
        'CLNREVSTAT': df.get('CLNREVSTAT', pd.Series(['.'] * len(df))),
        'CLNSIG': df.get('CLNSIG', pd.Series(['.'] * len(df))),
        'ONCDN': df.get('ONCDN', pd.Series(['.'] * len(df))),
        'ONCDISDB': df.get('ONCDISDB', pd.Series(['.'] * len(df))),
        'ONCREVSTAT': df.get('ONCREVSTAT', pd.Series(['.'] * len(df))),
        'ONC': df.get('ONC', pd.Series(['.'] * len(df))),
        'SCIDN': df.get('SCIDN', pd.Series(['.'] * len(df))),
        'SCIDISDB': df.get('SCIDISDB', pd.Series(['.'] * len(df))),
        'SCIREVSTAT': df.get('SCIREVSTAT', pd.Series(['.'] * len(df))),
        'SCI': df.get('SCI', pd.Series(['.'] * len(df)))
    }
    
    output_columns1 = {
        'Chr': df['Chr'],
        'Start': df['Start'],
        'End': df['End'],
        'Ref': df['Ref'],
        'Alt': df['Alt'],
        'AF': df.apply(lambda row: extract_vcf_info(row, gt_column), axis=1)['AF'],
        'DP': df.apply(lambda row: extract_vcf_info(row, gt_column), axis=1)['DP'],
        'Depth_Ref': df.apply(lambda row: extract_vcf_info(row, gt_column), axis=1)['Depth_Ref'],
        'Depth_Alt': df.apply(lambda row: extract_vcf_info(row, gt_column), axis=1)['Depth_Alt'],
        'Genotype': df.apply(lambda row: extract_vcf_info(row, gt_column), axis=1)['Genotype']
    }

    if gt2_column is not None:
        output_columns1.update({
            'AF_2': df.apply(lambda row: extract_vcf_info(row, gt2_column), axis=1)['AF'],
            'DP_2': df.apply(lambda row: extract_vcf_info(row, gt2_column), axis=1)['DP'],
            'Depth_Ref_2': df.apply(lambda row: extract_vcf_info(row, gt2_column), axis=1)['Depth_Ref'],
            'Depth_Alt_2': df.apply(lambda row: extract_vcf_info(row, gt2_column), axis=1)['Depth_Alt'],
            'Genotype_2': df.apply(lambda row: extract_vcf_info(row, gt2_column), axis=1)['Genotype']
        })

    output_columns1.update({
        'rsID': df.get('Otherinfo6', pd.Series(['.'] * len(df))),
        'AF_exome': df.apply(extract_gnomad_fields, axis=1)['AF_exome'],
        'AF_popmax_exome': df.apply(extract_gnomad_fields, axis=1)['AF_popmax_exome'],
        'AF_genome': df.apply(extract_gnomad_fields, axis=1)['AF_genome'],
        'AF_popmax_genome': df.apply(extract_gnomad_fields, axis=1)['AF_popmax_genome'],
        'TIER': df.apply(extract_gnomad_fields, axis=1)['TIER']
    })

    # Add optional columns
    optional_columns = ['pHaplo', 'pTriplo', '_loeuf ', '_pli']
    optional_output_names = ['pHap', 'pTrio', 'leouf', 'pLI']
    for input_col, output_col in zip(optional_columns, optional_output_names):
        if input_col in df.columns:
            output_columns1[output_col] = df[input_col]
        else:
            output_columns1[output_col] = pd.Series(['.'] * len(df))

    # Create output DataFrame
    output_df = pd.DataFrame({**output_columns1, **output_columns2})

    # Append remaining columns
    used_columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 
                    'CLNREVSTAT', 'CLNSIG', 'ONCDN', 'ONCDISDB', 'ONCREVSTAT', 'ONC', 
                    'SCIDN', 'SCIDISDB', 'SCIREVSTAT', 'SCI', 'Tier'] + optional_columns
    excluded_columns = used_columns + ['geneName', 'Gene']
    remaining_columns = [col for col in df.columns if col not in excluded_columns]

    if remaining_columns:
        remaining_df = df[remaining_columns]
        output_df = pd.concat([output_df, remaining_df], axis=1)

    # Write to output file
    output_df.to_csv(output_file, sep="\t", index=False)
    print(f"Reformatted data written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reformat a VCF-like dataset into a tab-delimited file with specified columns."
    )
    parser.add_argument("input_file", type=str, help="Path to the input tab-delimited file")
    parser.add_argument("output_file", type=str, help="Path to the output tab-delimited file")
    parser.add_argument("--info-column", type=str, default="Otherinfo12", 
                        help="Name of the column containing VCF field order (default: Otherinfo12)")
    parser.add_argument("--gt-column", type=str, default="Otherinfo13", 
                        help="Name of the column containing VCF genotype field (default: Otherinfo13)")
    parser.add_argument("--gt2-column", type=str, default=None, 
                        help="Name of the column containing second VCF genotype field (optional, default: None)")
    parser.add_argument("--gnomad-column", type=str, default="Otherinfo11", 
                        help="Name of the column containing VCF INFO fields (default: Otherinfo11)")
    
    args = parser.parse_args()
    reformat_vcf_dataset(args.input_file, args.output_file, args.info_column, args.gt_column, 
                         args.gt2_column, args.gnomad_column)
