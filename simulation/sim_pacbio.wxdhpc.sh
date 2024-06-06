#!/bin/bash
#SBATCH --partition=tiny
#SBATCH --job-name="batch_sim"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --ntasks-per-node=1
#SBATCH --output=log.batchsim.log
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00

set -e

# Directory for output
outdir=batch_sim_pacbio_CLR
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

eval "$(conda shell.bash hook)"
conda activate speclong
# Error rates from 0.78 to 0.9 with step 0.01
for err in $(seq 0.78 0.01 0.9); do
    # Depth from 5 to 50 with step 5
    for depth in $(seq 5 5 50); do
        sample=pacbio_dp${depth}_acc${err}

        # Create sample directory
        mkdir -p $outdir/$sample
        echo "gen txt file for $sample"

        # Simulate data
        perl simu.data.haplotype.hla.pl $sample $outdir/$sample
        echo """
        pbsim --prefix $outdir/$sample/$sample \
              --strategy wgs \
              --genome $outdir/$sample/$sample.HLA.sep.fa \
              --depth $depth \
              --method qshmm \
              --qshmm pbsim3_model/QSHMM-RSII.model \
              --accuracy-mean $err

        """
        # Run PBSIM2 to generate CLR reads
        pbsim --prefix $outdir/$sample/$sample \
              --strategy wgs \
              --genome $outdir/$sample/$sample.HLA.sep.fa \
              --depth $depth \
              --method qshmm \
              --qshmm pbsim3_model/QSHMM-RSII.model \
              --accuracy-mean $err

        # Combine and compress fastq files
        cat $outdir/$sample/${sample}_*.fastq > $outdir/$sample/${sample}.fastq
        rm $outdir/$sample/${sample}_*.fastq
        rm $outdir/$sample/${sample}_*.ref
        rm $outdir/$sample/${sample}_*.maf
        gzip -f $outdir/$sample/${sample}.fastq

        # Run main analysis script
        python3 ../scripts/main.py -n $sample -o $outdir -j 15 -y pacbio -i HLA -r $outdir/$sample/${sample}.fastq.gz --db ../db/ --min_identity $err

        # Evaluation
        python3 ../evaluation/assess_typing.py -i HLA --true $outdir/$sample/$sample.HLA.hap.alleles.txt --infer $outdir/$sample/${sample}.HLA.type.result.txt > $outdir/$sample/acc.step1.log
        python3 ../evaluation/assess_typing.py -i HLA --true $outdir/$sample/$sample.HLA.hap.alleles.txt --infer $outdir/$sample/hlala.like.results.txt > $outdir/$sample/acc.step2.log
    done
done