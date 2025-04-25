## given a fasta file., change the name of the contig
## using biopython

from Bio import SeqIO
import sys
import os

def change_db_name(fasta_file):
    """
    Change the name of the contig in a fasta file to a new name.
    """
    # Check if the file exists
    if not os.path.isfile(fasta_file):
        print(f"File {fasta_file} does not exist.")
        return

    # Read the fasta file
    records = SeqIO.parse(fasta_file, "fasta")
    
    # Create a new list to store the modified records
    modified_records = []
    
    # Iterate through each record and change the name
    for record in records:
        field = record.description.split()
        new_name = field[1]
        record.id = new_name
        modified_records.append(record)
    
    # Write the modified records to a new fasta file
    output_file = os.path.splitext(fasta_file)[0] + "_modified.fasta"
    SeqIO.write(modified_records, output_file, "fasta")
    
    print(f"Modified fasta file saved as {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python change_db_name.py <fasta_file> ")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    
    change_db_name(fasta_file)