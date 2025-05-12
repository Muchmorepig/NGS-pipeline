#!/usr/bin/env python3
"""
Script to extract protein sequences from a GenBank file and output them in FASTA format
Usage: python gbff2fasta.py input.gb output_file
"""

import sys
import re

def extract_proteins_to_fasta(input_file, output_file):
    """
    Extract protein sequences from GenBank file and save to a single FASTA file
    """
    # Initialize variables
    protein_sequences = []
    current_gene = None
    current_protein = None
    current_product = None
    translation_text = ""
    
    # First, read the entire file content
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Find all CDS sections with translations
    cds_pattern = re.compile(r'\s+CDS\s+.*?\n(.*?)/translation="(.*?)"', re.DOTALL)
    cds_matches = cds_pattern.findall(content)
    
    for cds_info, translation in cds_matches:
        # Extract gene name
        gene_match = re.search(r'/gene="(.*?)"', cds_info)
        if gene_match:
            current_gene = gene_match.group(1)
        else:
            continue  # Skip if no gene name found
        
        # Extract protein ID
        protein_match = re.search(r'/protein_id="(.*?)"', cds_info)
        if protein_match:
            current_protein = protein_match.group(1)
        else:
            continue  # Skip if no protein ID found
        
        # Extract product description
        product_match = re.search(r'/product="(.*?)"', cds_info)
        if product_match:
            current_product = product_match.group(1)
        else:
            current_product = None
        
        # Clean up translation text (remove spaces and newlines)
        clean_translation = re.sub(r'\s+', '', translation)
        
        # Store the protein sequence
        protein_sequences.append({
            'gene': current_gene,
            'protein_id': current_protein,
            'product': current_product,
            'sequence': clean_translation
        })
    
    # Write all protein sequences to a single FASTA file
    with open(output_file, 'w') as f_out:
        for protein in protein_sequences:
            header = f">{protein['gene']} protein_id={protein['protein_id']}"
            if protein['product']:
                header += f" {protein['product']}"
            f_out.write(f"{header}\n")
            
            # Write sequence in lines of 60 characters
            sequence = protein['sequence']
            for i in range(0, len(sequence), 60):
                f_out.write(f"{sequence[i:i+60]}\n")
    
    print(f"Extracted {len(protein_sequences)} protein sequences to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gbff2fasta.py input.gb output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    extract_proteins_to_fasta(input_file, output_file)