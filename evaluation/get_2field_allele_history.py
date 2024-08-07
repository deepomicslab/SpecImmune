import re

def convert_to_two_field(allele):
    """
    Convert a high-field allele notation to a two-field allele notation.
    """
    match = re.match(r'^(\w+\*\d+:\d+)', allele)
    if match:
        return match.group(1)
    return allele

def process_allele_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or line.strip() == "":
                outfile.write(line)
            else:
                parts = line.strip().split(',')
                hla_id = parts[0]
                alleles = parts[1:]
                two_field_alleles = [convert_to_two_field(allele) for allele in alleles]
                outfile.write(f"{hla_id},{','.join(two_field_alleles)}\n")

input_file = 'Allelelist_history.txt'
output_file = 'Allelelist_history_2field.txt'
process_allele_file(input_file, output_file)