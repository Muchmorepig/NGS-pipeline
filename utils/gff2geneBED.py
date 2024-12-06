#!/usr/bin/env python3

import sys
import re

def parse_attributes(attr_string):
    """Parse the attributes column of GFF file."""
    attrs = {}
    for attribute in attr_string.strip(';').split(';'):
        if '=' in attribute:
            key, value = attribute.split('=')
            attrs[key] = value
    return attrs

def gff_to_bed(gff_file, bed_file):
    """Convert GFF file to BED format, extracting only gene features."""
    with open(gff_file, 'r') as infile, open(bed_file, 'w') as outfile:
        for line in infile:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Parse GFF fields
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
                
            chromosome = fields[0]
            feature_type = fields[2]
            start = int(fields[3]) - 1  # Convert to 0-based coordinate system
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_attributes(fields[8])
            
            # Only process gene features
            if feature_type == 'gene':
                gene_id = attributes.get('ID', 'NA')
                gene_name = attributes.get('Name', gene_id)
                
                # Write in BED format: chr start end name score strand
                bed_line = f"{chromosome}\t{start}\t{end}\t{gene_name}\t0\t{strand}\n"
                outfile.write(bed_line)

if __name__ == "__main__":
    # Check if input and output filenames are provided
    if len(sys.argv) != 3:
        print("Usage: python gff2bed.py input.gff output.bed")
        sys.exit(1)
    
    input_gff = sys.argv[1]
    output_bed = sys.argv[2]
    
    # Convert GFF to BED
    gff_to_bed(input_gff, output_bed)
    print(f"Conversion complete. BED file saved as: {output_bed}")
