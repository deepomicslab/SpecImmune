# set -e 
## use pbsim2


outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads3/
resultdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results3/

for dp in 5 10 15 20 25 30 
do
for i in {21..50}
do
    sample=KIR_dp${dp}_acc90_$i
    echo $sample

    #### simulation
    mkdir $outdir/$sample
    perl simu.data.haplotype.kir2.pl $sample $outdir/$sample
    pbsim --prefix $outdir/$sample/$sample --depth $dp --hmm_model pbsim_model/P4C2.model --accuracy-mean 0.9 $outdir/$sample/$sample.KIR.sep.fa

    cat $outdir/$sample/${sample}_*fastq>$outdir/$sample/${sample}.fastq
    rm $outdir/$sample/${sample}_*fastq
    rm $outdir/$sample/${sample}_*.ref
    rm $outdir/$sample/${sample}_*.maf
    gzip -f $outdir/$sample/${sample}.fastq


    #### run
    python3 ../scripts/main.py -n $sample -o $resultdir -j 15 -y pacbio -i KIR -r $outdir/$sample/${sample}.fastq.gz -k 1 --hete_p 0.2

    # #### evaluation
    # python3 ../evaluation/assess_read_bin.py $resultdir/$sample/${sample}.read_binning.txt $outdir/$sample/$sample.KIR.sep.fa $outdir/$sample/${sample}.fastq.gz 
    # python3 ../evaluation/assess_typing.py -i KIR --true $outdir/$sample/$sample.KIR.hap.alleles.txt --infer $resultdir/$sample/${sample}.KIR.type.result.txt 
    # echo "----split-----$sample"
    # python3 ../evaluation/assess_typing.py -i KIR --true $outdir/$sample/$sample.KIR.hap.alleles.txt --infer $resultdir/$sample/${sample}.KIR.final.type.result.txt 
    # break

done
done

# outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads3/
# resultdir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results3/
# sample=KIR_dp40_acc90_14
# # python3 ../scripts/main.py -n $sample -o $resultdir -j 15 -y pacbio -i KIR -r $outdir/$sample/${sample}.fastq.gz -k 1 --hete_p 0.2
    
#     python3 ../evaluation/assess_read_bin.py $resultdir/$sample/${sample}.read_binning.txt $outdir/$sample/$sample.KIR.sep.fa $outdir/$sample/${sample}.fastq.gz 
#     # python3 ../evaluation/assess_typing.py -i KIR --true $outdir/$sample/$sample.KIR.hap.alleles.txt --infer $resultdir/$sample/${sample}.KIR.type.result.txt 
#     echo "----split-----$sample"
#     python3 ../evaluation/assess_typing.py -i KIR --true $outdir/$sample/$sample.KIR.hap.alleles.txt --infer $resultdir/$sample/${sample}.KIR.final.type.result.txt 

#     # result=$resultdir/$sample/$sample.KIR.final.type.result.txt
#     # ## check the line number of result, if the line number is 2, then the result is empty
#     # if [ -f $result ]; then
#     #     line_num=$(wc -l $result | awk '{print $1}')
#     #     if [ $line_num -eq 2 ]; then
#     #         echo "handle $sample"
#     #         python3 ../scripts/main.py -n $sample -o $resultdir -j 15 -y pacbio -i KIR -r $outdir/$sample/${sample}.fastq.gz -k 1 --hete_p 0.2
#     #         # continue
#     #     fi
#     # fi
