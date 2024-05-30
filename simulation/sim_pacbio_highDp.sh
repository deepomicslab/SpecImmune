for i in {1..2}
do
sample=child_$i
mkdir high_pacbio/$sample
pbsim --data-type CLR --seed 88 --accuracy-mean 0.85 --accuracy-min 0.80 --prefix high_pacbio/$sample/$sample --depth 200 --model_qc pacbio/model_qc_clr fasta/$sample.fasta
cat high_pacbio/$sample/${sample}_*fastq>high_pacbio/$sample/${sample}.fastq
rm high_pacbio/$sample/${sample}_*fastq
done
