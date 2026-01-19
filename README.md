# Pipeline for generation of synonymous mutant library

## Introduction

This project aims to systematically generate a library of synonymous mutants in the ddlB gene, a protein model known to undergo entanglement during its folding process. Synonymous codons, although encoding the same amino acid, can influence the folding trajectory of proteins. However, pinpointing which codons or regions cause these folding changes remains a significant challenge.

The ddlB gene encodes a protein of 306 amino acids. The mutagenesis focuses on the region spanning amino acids 97 (lysine) to 231 (alanine), covering 135 codons (405 nucleotides). Notably, the entanglement event occurs between amino acids 98 (leucine) and 146 (leucine), making this region critical for studying folding dynamics.

The goal is to create a comprehensive library of synonymous variants using oligopools. Each oligopool targets a stretch of codons scanned with single synonymous mutations. Additionally, random second- and third-order combinations of synonymous mutations will be generated for each stretch to explore multi-mutant effects.
A key aspect of this project is to track the cloned libraries throughout the experimental pipeline. Automated liquid handling will be employed, necessitating rigorous sample and variant tracking to ensure data integrity and reproducibility.

## Files in the folder

degenerated_codons.xlsx: Contains a dataframe of degenerated codons per amino acid. For example, Alanine (Ala) is encoded by GCN.

## Pipeline

1. Generation of  oligos: In each position of an oligo, the other mutant will be tested for every position. For that, uses N, Y or R to meke them degenerated.

