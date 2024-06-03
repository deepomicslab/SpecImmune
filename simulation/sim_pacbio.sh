## use pbsim2


outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads/


for i in {1..5}
do
    sample=pacbio_dp50_acc98_$i

    #### simulation
    # mkdir $outdir/$sample
    # perl simu.data.haplotype.hla.pl $sample $outdir/$sample
    # pbsim --prefix $outdir/$sample/$sample --depth 50 --hmm_model pbsim_model/P4C2.model --accuracy-mean 0.98 $outdir/$sample/$sample.HLA.sep.fa

    # cat $outdir/$sample/${sample}_*fastq>$outdir/$sample/${sample}.fastq
    # rm $outdir/$sample/${sample}_*fastq
    # rm $outdir/$sample/${sample}_*.ref
    # rm $outdir/$sample/${sample}_*.maf
    # gzip -f $outdir/$sample/${sample}.fastq


    #### run
    python ../scripts/main.py -n $sample -o $outdir -j 15 -y pacbio -i HLA -r $outdir/$sample/${sample}.fastq.gz --db ../db/ 

    #### evaluation
    break

done



