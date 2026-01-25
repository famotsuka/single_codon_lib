from scripts.design_lib import gene, generate_fragments, single_mutants, oligo_pre_sequences, mutate_sequence, cut_fragment

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

    nb_of_frags = 3 # number of fragments to divide the gene into for the golden gate fragment

    frag_5prime_bsaI_site = "ACTGTAGGTCTCCT"  # BsaI site for 5'

    '''
    Core of the script
    '''
    # 1) Load gene information and create the codon dataframe with single mutants per position
    gene_seq, gene_id, codon_df = gene(gene_of_interest) #apply the function to extract the gene information
    codon_df = single_mutants(codon_df)  # Generate single mutants and add to dataframe
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

        # Save temp_df without mut_gene column
        output_df = temp_df.drop(columns=['mut_gene'])
        output_file = f"./outputs/{frag_name}_oligos.csv"
        output_df.to_csv(output_file, index=False)
        print(f"\nSaved {frag_name} to {output_file}")



if __name__ == "__main__":
    main()
