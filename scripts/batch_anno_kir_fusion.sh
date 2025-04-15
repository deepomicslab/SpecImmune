

samples="HG01489 HG01915 HG03074 HG03097 HG03133 HG03247 HG03470 HG03484 HG03485 HG03511 HG03539 NA18873 NA19096 NA19327 NA19328 NA19347"

for sample in $samples
do
    echo $sample
    python kir_annotation.py -i speclong_kir_out/${sample}/${sample}/${sample}.KIR.final.type.result.formatted.txt -o speclong_kir_out/${sample}/${sample}/${sample}.KIR.final.type.result.formatted.anno.txt
done