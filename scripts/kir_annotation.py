import pandas as pd
import argparse
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description="Annotate KIR3DL1-KIR3DL2 fusion in KIR typing results.")
    parser.add_argument("-i", "--input", required=True, help="Input KIR typing TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output annotated TSV file")
    return parser.parse_args()

def extract_prefix(allele):
    """Extract prefix like KIR3DL1*059 from full allele like KIR3DL1*05901"""
    match = re.match(r"(KIR3DL1\*\d{3})", allele)
    return match.group(1) if match else None

def main():
    args = parse_arguments()
    input_file = args.input
    output_file = args.output

    # Define fusion allele prefixes
    fusion_prefixes = {"KIR3DL1*059", "KIR3DL1*060", "KIR3DL1*061"}

    # Read file and skip metadata lines
    with open(input_file) as f:
        lines = f.readlines()

    header_line_index = next(i for i, line in enumerate(lines) if line.startswith("Locus"))
    header = lines[header_line_index].strip().split('\t')
    data_lines = lines[header_line_index + 1:]
    data = [line.strip().split('\t') for line in data_lines]

    df = pd.DataFrame(data, columns=header)
    df["Fusion"] = ""

    # Track fusion alleles by chromosome index (e.g., '1' or '2')
    fusion_allele_map = {}

    # Find KIR3DL1 fusion alleles
    kir3dl1_df = df[df["Locus"] == "KIR3DL1"]
    for _, row in kir3dl1_df.iterrows():
        chr_index = row["Chromosome"]
        alleles = row["Genotype"].split(";")
        for allele in alleles:
            prefix = extract_prefix(allele.strip())
            if prefix in fusion_prefixes:
                fusion_allele_map.setdefault(chr_index, []).append(allele.strip())

    # Annotate KIR3DL2 rows
    kir3dl2_indices = df[df["Locus"] == "KIR3DL2"].index
    for idx in kir3dl2_indices:
        chr_index = df.at[idx, "Chromosome"]
        if chr_index in fusion_allele_map:
            fusion_alleles = ";".join(fusion_allele_map[chr_index])
            df.at[idx, "Fusion"] = fusion_alleles
            df.at[idx, "One_guess"] = fusion_alleles
            df.at[idx, "Genotype"] = fusion_alleles  # Optional: comment this line out if not needed

    # Save annotated file
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Fusion annotation complete. Output written to: {output_file}")

if __name__ == "__main__":
    main()