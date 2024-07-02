set -e 
## use pbsim2


outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads/
resultdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/


for i in {1..5}
do
    sample=IG_TR_dp50_acc98_$i

    #### simulation
    # python simu.data.haplotype.VDJ.py $sample $outdir/$sample

    # pbsim --prefix $outdir/$sample/$sample --depth 50 --hmm_model pbsim_model/P4C2.model --accuracy-mean 0.95 $outdir/$sample/$sample.IG.TR.hap1.fa
    # pbsim --prefix $outdir/$sample/$sample --depth 50 --hmm_model pbsim_model/P4C2.model --accuracy-mean 0.95 $outdir/$sample/$sample.IG.TR.hap2.fa
    # cat $outdir/$sample/${sample}_*fastq>$outdir/$sample/${sample}.fastq
    # rm $outdir/$sample/${sample}_*fastq
    # rm $outdir/$sample/${sample}_*.ref
    # rm $outdir/$sample/${sample}_*.maf
    # gzip -f $outdir/$sample/${sample}.fastq


    #### run
    python3 ../scripts/main.py -n $sample -o $resultdir -j 15 -y pacbio -i IG_TR -r $outdir/$sample/${sample}.fastq.gz 

    #### evaluation
    # python3 ../evaluation/assess_typing.py -i IG_TR --true $outdir/$sample/$sample.IG_TR.hap.alleles.txt --infer $outdir/$sample/${sample}.HLA.type.result.txt 
    break

done



