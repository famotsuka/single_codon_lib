import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq  

'''
Description
Contains the functions for generating a library of oligopools 
for single synonymous site saturation.
Use "main.py"
'''

def gene(gene_fasta):
    '''
    Read gene sequence from fasta file and create a dataframe with codon information.
    Returns the gene sequence and a dataframe with codon position, codon, and amino acid.
    '''    
    for record in SeqIO.parse(gene_fasta, "fasta"):
        gene_seq = str(record.seq)
        gene_id = record.id
    
    # Create the WT codon dataframe
    codons_data = []
    for i in range(0, len(gene_seq), 3):
        codon = gene_seq[i:i+3]
        if len(codon) == 3:  # Only include complete codons
            seq_obj = Seq(codon)
            aa = seq_obj.translate()
            position = (i // 3) + 1  # Position starting from 1
            codons_data.append({
                'codon_position': position,
                'wt_codon': codon,
                'amino_acid': str(aa)
            })
    
    codon_df = pd.DataFrame(codons_data)
    return gene_seq, gene_id, codon_df

def generate_fragments(data_frame, nb_of_frags, start_window, end_window):
    '''
    Divide the gene into a specified number of fragments for independent golden gate assembly.
    Returns a dictionary with fragment keys (frag1, frag2, etc.) containing codon positions.
    '''
    fragments = {}
    # Filter rows within the window
    window_df = data_frame[(data_frame['codon_position'] >= start_window) & 
                           (data_frame['codon_position'] <= end_window)]
    codons_per_frag = len(window_df) // nb_of_frags
    print(f"\nNumber of codons per fragment: {codons_per_frag} \n")

    # Divide window_df into fragments
    for i in range(nb_of_frags):
        start_index = i * codons_per_frag
        if i == nb_of_frags - 1:  # Last fragment takes the remainder
            end_index = len(window_df)
        else:
            end_index = (i + 1) * codons_per_frag
        
        fragment_data = window_df.iloc[start_index:end_index]
        frag_name = f"frag{i+1}"
        fragments[frag_name] = fragment_data['codon_position'].tolist()

    return fragments

def single_mutants(codon_df, codon_to_rarity):
    '''
    Generate a list of synonymous codons for each amino acid (excluding the WT codon).
    Adds a new column called 'single_muts' to the dataframe.
    Adds a new column called 'single_muts_rarity' to the dataframe to track codon mutant rarity based on the codon usage table.
    Returns the modified dataframe.
    '''
    from Bio.Data import CodonTable
    
    # Get the standard codon table
    codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
    
    # amino acid -> list of codons
    aa_to_codons = {}
    for codon, aa in codon_table.forward_table.items():
        aa_to_codons.setdefault(aa, []).append(codon)
    aa_to_codons['*'] = list(codon_table.stop_codons)

    def get_synonymous_and_rarity(aa, wt_codon):
        if aa not in aa_to_codons:
            return [], []
        muts = [c for c in aa_to_codons[aa] if c != wt_codon]
        rarities = [codon_to_rarity.get(c, None) for c in muts]
        return muts, rarities

    codon_df[["single_muts", "single_muts_rarity"]] = codon_df.apply(
        lambda row: pd.Series(
            get_synonymous_and_rarity(row["amino_acid"], row["wt_codon"])
        ),
        axis=1,
    )
    
    return codon_df

def mutate_sequence(wt_seq, codon_pos, new_codon):
    '''
    Generate the full mutated sequence of the gene.
    '''
    start = (codon_pos - 1) * 3
    end = start + 3
    return wt_seq[:start] + new_codon + wt_seq[end:]

def cut_fragment(seq, start_pos, end_pos):
    '''
    Cut a fragment from the sequence based on start and end codon positions of the fragment
    and extend it 2 nucleotides before and 9 nucleotides after.
    '''
    start = (start_pos - 1) * 3 - 6  # Extend 6 nucleotides before to give space for BsaI overhang.
    end = (end_pos - 1) * 3 + 3 + 12  # Extend 12 nucleotides after for ensure annealing at 3' end.
    
    # Ensure we don't go out of bounds
    start = max(0, start)
    end = min(len(seq), end)
    
    return seq[start:end]

def oligo_pre_sequences(gene_seq, temp_df):
    '''
    Generate the pre-sequence for oligos using the single_mutants.
    For each row, creates mutant sequences by replacing the wt-codon with each single_mut.
    Cuts fragments from each mutant sequence using cut_fragment.
    Returns the dataframe with 'mut_gene' and 'mut_frag' columns.
    '''
    mut_gene_list = []
    mut_frag_list = []
    
    # Get start and end positions for the fragment
    start_pos = temp_df['codon_position'].iloc[0]
    end_pos = temp_df['codon_position'].iloc[-1]
    
    for _, row in temp_df.iterrows():
        codon_pos = row["codon_position"]
        mut_sequences = []
        mut_fragments = []
        
        for mut in row["single_muts"]:
            mutated_seq = mutate_sequence(gene_seq, codon_pos, mut)
            mut_sequences.append(mutated_seq)
            
            # Cut fragment from the mutated sequence
            fragment = cut_fragment(mutated_seq, start_pos, end_pos)
            mut_fragments.append(fragment)
        
        mut_gene_list.append(mut_sequences)
        mut_frag_list.append(mut_fragments)
    
    temp_df['mut_gene'] = mut_gene_list
    temp_df['mut_frag'] = mut_frag_list
    return temp_df

def split_rarity(df):
    """
    Split a dataframe with list-columns into two:
    - df_rare: keeps only entries where rarity == 'rare'
    - df_common: keeps only entries where rarity == 'common'
    
    Applies filtering consistently across:
        - single_muts
        - single_muts_rarity
        - mut_frag
    """
    rare_rows = []
    common_rows = []

    for _, row in df.iterrows():
        muts = row["single_muts"]
        rarities = row["single_muts_rarity"]
        frags = row["mut_frag"]

        # Build index masks
        rare_idx = [i for i, r in enumerate(rarities) if r == "rare"]
        common_idx = [i for i, r in enumerate(rarities) if r == "common"]

        # Build filtered rows
        if rare_idx:
            rare_rows.append({
                **row.to_dict(),
                "single_muts": [muts[i] for i in rare_idx],
                "single_muts_rarity": [rarities[i] for i in rare_idx],
                "mut_frag": [frags[i] for i in rare_idx],
            })

        if common_idx:
            common_rows.append({
                **row.to_dict(),
                "single_muts": [muts[i] for i in common_idx],
                "single_muts_rarity": [rarities[i] for i in common_idx],
                "mut_frag": [frags[i] for i in common_idx],
            })

    df_rare = pd.DataFrame(rare_rows)
    df_common = pd.DataFrame(common_rows)

    return df_rare, df_common

def explode_unique_sequences(df):
    # Explode mut_frag into individual rows
    exploded = df.explode('mut_frag')
    
    # Drop duplicates to keep only unique (Pool name, sequence) pairs
    exploded = exploded[['Pool name', 'mut_frag']].drop_duplicates()
    
    # Rename mut_frag to sequence
    exploded.rename(columns={'mut_frag': 'sequence'}, inplace=True)
    
    return exploded.reset_index(drop=True)


if __name__ == "__main__":
    main()
