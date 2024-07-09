from Bio import SeqIO
import gzip

# Define the path to your FASTQ file and the read name you're looking for
fastq_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/SRR19372746.fastq.gz"
read_name = "SRR19372746.149369"

# Open the gzipped FASTQ file and search for the read
found = False
with gzip.open(fastq_file, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        if record.id == read_name:
            print(f"Read name: {record.id}")
            print(f"Sequence: {record.seq}")
            # print(f"Quality: {record.letter_annotations['phred_quality']}")
            found = True
            break

# if not found:
#     print("Read name not found in the FASTQ file.")
#     if record.id == read_name:
#         print(f"Read name: {record.id}")
#         print(f"Sequence: {record.seq}")
#         print(f"Quality: {record.letter_annotations['phred_quality']}")
#         found = True
        # break

