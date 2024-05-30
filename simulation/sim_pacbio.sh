## use pbsim2


outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads/


for i in {1..5}
do
    sample=pacbio_dp50_acc98_$i
    mkdir $outdir/$sample
    # pbsim --data-type CLR --seed 88 --accuracy-mean 0.85 --accuracy-min 0.80 --prefix $outdir/$sample/$sample --depth 10 --model_qc $outdir/model_qc_clr fasta/$sample.fasta
    perl simu.data.haplotype.hla.pl $sample $outdir/$sample
    pbsim --prefix $outdir/$sample/$sample --depth 50 --hmm_model pbsim_model/P4C2.model --accuracy-mean 0.98 $outdir/$sample/$sample.HLA.sep.fa

    cat $outdir/$sample/${sample}_*fastq>$outdir/$sample/${sample}.fastq
    rm $outdir/$sample/${sample}_*fastq
    rm $outdir/$sample/${sample}_*.ref
    rm $outdir/$sample/${sample}_*.maf
    gzip -f $outdir/$sample/${sample}.fastq
    # break
done



