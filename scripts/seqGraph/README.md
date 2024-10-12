# seqGraph
SEQ is a conjugated assembly algorithm that constructs a conjugate map according to Barcode connection information between Contig, and finds the most likely connection relationship between Contig.
## Requirement
htslib
## Install
```
git clone https://github.com/panguangze/seqGraph
mkdir build
cd build
cmake .. && make
```
## Running
#### make input
```
Add tgs data into assembly process
p_tgs -b tgs.bam -c contig.fasta -o out.fasta
```
tgs.bam aligned the tgs reads to the assemblied raw contigs.
```
overlap_end input.bam graph.txt 500,
```
The input.bam is aligned NGS reads to the contigs. The last parameters is the end length of contig for barcode statistic
```
seqGraph \
    -g  graph.txt\ // graph input
    -o  assembly_orders.txt\ //assemly contigs order and orientations
    -i  1 //iteration times, more iteration times longer scaf length but less accuracy
```
make fasta from assembly result
```
python scripts/make_fa.py assembly_orders.txt out.fasta
```
