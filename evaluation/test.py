import csv
from collections import defaultdict

def read_tsv_to_dict(file_path):
    # Initialize the dictionary to hold the parsed data
    data_dict = defaultdict(dict)
    
    with open(file_path, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        for row in reader:
            sample = row['Sample']
            # For each gene, add it to the sample's dictionary
            for gene in row:
                if gene != 'Sample':
                    data_dict[sample][gene] = row[gene]
    
    return data_dict

# Example usage
file_path = 'HLA-LA.merge.result.txt'  # Replace with the path to your actual TSV file
data_dict = read_tsv_to_dict(file_path)

# Print the resulting dictionary
import pprint
pprint.pprint(data_dict)