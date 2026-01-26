# Pipeline for generation of synonymous mutant library

## Introduction

This project aims to systematically generate a library of synonymous mutants in the ddlB gene, a protein model known to undergo entanglement during its folding process. Synonymous codons, although encoding the same amino acid, can influence the folding trajectory of proteins. However, pinpointing which codons or regions cause these folding changes remains a significant challenge.

The ddlB gene encodes a protein of 306 amino acids. The mutagenesis focuses on the region spanning amino acids 97 (lysine) to 234 (glutamate), covering 138 codons (414 nucleotides). The entanglement event occurs between amino acids 98 (leucine) and 146 (leucine), making this region critical for studying folding dynamics.

The goal is to create a library of synonymous variants using oligopools. Each oligopool targets a stretch of codons scanned with single synonymous mutations. Additionally, random second- and third-order combinations of synonymous mutations might be addded later to explore multi-mutant or genetic background effects.

## Files in the repository

- "degenerated_codons.xlsx": Contains a dataframe of degenerated codons per amino acid. For example, Alanine (Ala) is encoded by GCN.

- "main.py": Is the main script to generate the site-saturation single synonymous library in a given gene.

- "scripts/design_lib.py": contains all the functions by which "main.py" is calling 

- "scripts/get_fasta.py": Script to run alone and get all the sequences in a fasta format of fragments to order.

## Usage

Clone the repository in your local computer:
`git clone https://github.com/famotsuka/single_codon_lib.git`

Synchronize environment with uv:
`uv sync`

activate the environment:
`source .venv/bin/activate`

Create two folders:
`mkdir inputs`
`mkdir outputs`

Add the gene of interest in "inputs/" using a fasta format file and starting with the +1 frame.

Define parameters and inputs inside "main.py", such as:
gene_of_interest = "./inputs/ddlB_gBlock.fasta"  # As fasta file and without stop condon
start_window_frag_Fw = 97 # aminoacid position
end_window_frag_Fw = 234 # aminoacid position
nb_of_frags = 3 # number of fragments to divide the gene into for the golden gate fragment
frag_5prime_bsaI_site = "ACTGTAGGTCTCCT"  # BsaI site for 5'

Run the script:
`python3 main.py`

## Pipeline

1. Generation of  oligos: Pooled single synonymous site saturation oligos per stretches of the gene.

2. Order tubes: The output is ready to sendo for oreder in IDT.