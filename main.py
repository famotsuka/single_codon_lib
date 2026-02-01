from scripts.design_lib import gene, generate_fragments, single_mutants, split_rarity
from scripts.design_lib import oligo_pre_sequences, mutate_sequence, cut_fragment, explode_unique_sequences
import pandas as pd

'''
Description
Run this script to generate a library of oligopools for single synonymous site saturation.
'''

def main():
    print("\nScript for generation of single-codon-lib \n")
    '''
    Inputs of script.
    '''
    gene_of_interest = "./inputs/ddlB_gBlock.fasta"  # As fasta file and without stop condon

    start_window_frag_Fw = 97 # aminoacid position
    end_window_frag_Fw = 234 # aminoacid position

    nb_of_frags = 2 # number of fragments to divide the gene into for the golden gate fragment

    frag_5prime_bsaI_site = "ACTGTAGGTCTCC"  # BsaI site for 5'

    codon_usage_path = "./codon_tables/cocoputs_ecoli_CAI.csv"  # codon usage table for E. coli

    '''
    Core of the script
    '''
    # 1) Load gene information and create the codon dataframe with single mutants per position
    gene_seq, gene_id, codon_df = gene(gene_of_interest) #apply the function to extract the gene information
    
    # Load codon usage table
    codon_table = pd.read_csv(codon_usage_path, sep=";")
    codon_table = codon_table[["codon", "Rarity"]]
    codon_to_rarity = dict(zip(codon_table["codon"], codon_table["Rarity"]))
    #print(f"\nCodon usage table head: \n{codon_table.head().to_string(index=False)}")
    print(f"\ndictionary codon to rarity head: \n{dict(list(codon_to_rarity.items())[:5])}")

    codon_df = single_mutants(codon_df, codon_to_rarity)  # Generate single mutants and add to dataframe
    print(f"Gene of interest ID: {gene_id} \n")
    print(f"Gene sequence: {gene_seq} \n")
    print(f"Gene info head: \n{codon_df.head().to_string(index=False)}")
    print(f"Gene info tail: \n{codon_df.tail().to_string(index=False)}")

   
    # 2) Generate fragments
    fragments = generate_fragments(codon_df, nb_of_frags, start_window_frag_Fw, end_window_frag_Fw)
    print(f"\nPositions to mutated in fragments: \n{fragments}")

    # 3) Generate pre-sequences containing the single mutants
    for frag_name, codon_positions in fragments.items():
        # Create temporary dataframe for this fragment
        temp_df = codon_df[codon_df['codon_position'].isin(codon_positions)].copy()
        
        # Add Pool name column
        first_pos = temp_df['codon_position'].iloc[0]
        last_pos = temp_df['codon_position'].iloc[-1]
        temp_df.loc[:, "Pool name"] = f"{frag_name.upper()}_{first_pos}_to_{last_pos}"
        
        # Generate mutant gene sequences
        temp_df = oligo_pre_sequences(gene_seq, temp_df)

        # Add BsaI site to 5' end of fragment
        temp_df["mut_frag"] = temp_df["mut_frag"].apply(lambda lst: [frag_5prime_bsaI_site + seq for seq in lst])

        # Save temp_df without mut_gene column
        output_info_df = temp_df.drop(columns=['mut_gene'])
        df_rare, df_common = split_rarity(output_info_df)
        print(f"\nRarity split for {frag_name}: \nPositions with rare codons: {len(df_rare)} \nPositions with common codons: {len(df_common)}")
        # -----------------------------
        # 1) Save FULL pool
        # -----------------------------
        full_info_file = f"./outputs/{frag_name}_FULL_oligos_info.csv"
        output_info_df.to_csv(full_info_file, index=False)
        print(f"\nSaved full pool for {frag_name} to {full_info_file}")

        # Generate oligos to order (full)
        df_to_order_full = explode_unique_sequences(output_info_df)
        df_to_order_full = df_to_order_full[df_to_order_full["sequence"].str.len() > 0]

        full_order_file = f"./outputs/{frag_name}_FULL_oligos_to_order.csv"
        df_to_order_full.to_csv(full_order_file, index=False)
        print(f"Saved full oligo pool to {full_order_file}")

        # -----------------------------
        # 2) Save RARE pool
        # -----------------------------
        rare_info_file = f"./outputs/{frag_name}_rare_oligos_info.csv"
        df_rare.to_csv(rare_info_file, index=False)
        print(f"\nSaved RARE pool for {frag_name} to {rare_info_file}")

        # Generate oligos to order (rare)
        df_to_order_rare = explode_unique_sequences(df_rare)
        df_to_order_rare = df_to_order_rare[df_to_order_rare["sequence"].str.len() > 0]

        rare_order_file = f"./outputs/{frag_name}_rare_oligos_to_order.csv"
        df_to_order_rare.to_csv(rare_order_file, index=False)
        print(f"Saved RARE oligo pool to {rare_order_file}")

        # -----------------------------
        # 3) Save COMMON pool
        # -----------------------------
        common_info_file = f"./outputs/{frag_name}_common_oligos_info.csv"
        df_common.to_csv(common_info_file, index=False)
        print(f"\nSaved COMMON pool for {frag_name} to {common_info_file}")

        # Generate oligos to order (common)
        df_to_order_common = explode_unique_sequences(df_common)
        df_to_order_common = df_to_order_common[df_to_order_common["sequence"].str.len() > 0]

        common_order_file = f"./outputs/{frag_name}_common_oligos_to_order.csv"
        df_to_order_common.to_csv(common_order_file, index=False)
        print(f"Saved COMMON oligo pool to {common_order_file}")

if __name__ == "__main__":
    main()
