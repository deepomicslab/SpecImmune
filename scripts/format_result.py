import csv

def format_hla_file(input_file, output_file):
    # Read the input file
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        lines = list(reader)

    # Extract the version line and the header
    version_line = lines[0][0]
    header = lines[1]
    
    # Extract data lines
    data_lines = lines[2:]

    # Process the data lines according to the specified rules
    formatted_lines = []
    for parts in data_lines:
        locus, chromosome, genotype, match_info, reads_num, step1_type, one_guess = parts
        
        if genotype == 'NA' and step1_type != '-':
            genotype = step1_type
        
        formatted_line = [locus, chromosome, genotype, match_info, reads_num, step1_type, one_guess]
        formatted_lines.append(formatted_line)

    # Write the output to a new file
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow([version_line])
        writer.writerow(header)
        writer.writerows(formatted_lines)

if __name__ == "__main__":
    import argparse

    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Format HLA result file")
    parser.add_argument("input_file", help="The input file containing HLA results")
    parser.add_argument("output_file", help="The output file to save formatted results")

    args = parser.parse_args()

    # Run the formatting function
    format_hla_file(args.input_file, args.output_file)