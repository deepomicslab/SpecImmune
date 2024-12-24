
# HLA
bash SpecComplex.sh -n $sample_name -1 $fq1 -2 $fq2 -o $output_dir -i HLA
# CYP
bash SpecComplex.v1.sh -n $sample_name -1 $fq1 -2 $fq2 -o $output_dir -i CYP -m 10 -j 5 -q 10 -r 0.1
# KIR
bash SpecComplex.sh -n $sample_name -1 $fq1 -2 $fq2 -o $output_dir -i KIR
# vdj
bash run.wxs.vdj.sh -n $sample_name -1 $fq1 -2 $fq2
