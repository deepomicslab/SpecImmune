#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <input_dir>"
    exit 1
fi

input_dir="$1"


if [ ! -d "$input_dir" ]; then
    echo "Input directory $input_dir does not exist."
    exit 1
fi

for fasta_file in "$input_dir"/whole/*.fasta; do
    if [[ -f "$fasta_file" ]]; then
        db_name="${fasta_file%.fasta}"
        echo "Creating BLAST database for $fasta_file..."
        makeblastdb -in "$fasta_file" -dbtype nucl -out "$db_name"
        
        echo "BLAST database created at: $db_name"
    else
        echo "No FASTA files found in the directory."
    fi
done