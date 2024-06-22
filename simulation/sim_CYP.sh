set -e 
## use pbsim2


outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads/
resultdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/


for i in {1..5}
do
    sample=CYP_dp50_acc98_$i

    #### simulation
    # mkdir $outdir/$sample
    # python simu.data.haplotype.cyp.py $sample $outdir/$sample
    # pbsim --prefix $outdir/$sample/$sample --depth 50 --hmm_model pbsim_model/P4C2.model --accuracy-mean 0.95 $outdir/$sample/$sample.CYP.sep.fa

    # cat $outdir/$sample/${sample}_*fastq>$outdir/$sample/${sample}.fastq
    # rm $outdir/$sample/${sample}_*fastq
    # rm $outdir/$sample/${sample}_*.ref
    # rm $outdir/$sample/${sample}_*.maf
    # gzip -f $outdir/$sample/${sample}.fastq


    #### run
    # python3 ../scripts/main.py -n $sample -o $resultdir -j 15 -y pacbio -i CYP -r $outdir/$sample/${sample}.fastq.gz --db ../db/ 

    #### evaluation
    python3 ../evaluation/assess_read_bin.py $resultdir/$sample/${sample}.assign.txt $outdir/$sample/$sample.CYP.sep.fa $outdir/$sample/${sample}.fastq.gz 
    python3 ../evaluation/assess_typing.py -i CYP --true $outdir/$sample/$sample.CYP.hap.alleles.txt --infer $resultdir/$sample/${sample}.CYP.type.result.txt 
    break

done



