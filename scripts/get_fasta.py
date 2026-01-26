import csv
import os

def csv_to_fasta(input_csv, output_fasta):
    with open(input_csv, 'r') as csv_file, open(output_fasta, 'w') as fasta_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            fasta_file.write(f">{row['Pool name']}\n{row['sequence']}\n")

if __name__ == "__main__":
    # Directory where THIS script lives
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Outputs directory relative to script
    output_dir = os.path.join(script_dir, "..", "outputs")

    # Ensure the directory exists
    if not os.path.isdir(output_dir):
        raise FileNotFoundError(f"Outputs directory not found: {output_dir}")

    # Process files
    for file_name in os.listdir(output_dir):
        if file_name.endswith("frag3_oligos_to_order.csv"):
            input_csv_path = os.path.join(output_dir, file_name)
            output_fasta_path = os.path.join(output_dir, file_name.replace(".csv", ".fasta"))
            csv_to_fasta(input_csv_path, output_fasta_path)

    print("FASTA files created.")
